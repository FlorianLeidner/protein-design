import json
import argparse
import traceback

from pathlib import Path

import esm
import esm.inverse_folding

import torch
import torch.nn.functional as F
import torch.multiprocessing as mp

import MDAnalysis as mda

import numpy as np




def load_backbone_coords(pdb_path: Path, masked_atoms: str = None) -> tuple[np.ndarray, list[list]]:
    """
    Load a PDB with MDAnalysis and return:
        coords  : (N, 3, 3) float32 array of [N, CA, C] per residue
        seq     : list of one-letter amino-acid codes (length N)

    Only the backbone heavy atoms N / CA / C are extracted.
    Missing atoms for a residue are filled with NaN so that the residue
    index still lines up with the sequence.
    """

    THREE_TO_ONE = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        # common non-standard mappings
        "HSD": "H", "HSE": "H", "HSP": "H",
        "CYX": "C", "GLH": "E", "ASH": "D",
    }

    u = mda.Universe(pdb_path)

    protein = u.select_atoms("protein")
    residues = protein.residues

    if masked_atoms is not None:
        masked_atoms = u.select_atoms(masked_atoms)
    else:
        masked_atoms = ()

    coords = []
    
    seq = [[], []]

    for res in residues:
        resname = res.resname.upper()
        aa = THREE_TO_ONE.get(resname, "X")

        row = np.full((3, 3), 0., dtype=np.float32)

        for k, atom_name in enumerate(["N", "CA", "C"]):
            
            sel = res.atoms.select_atoms(f"name {atom_name}")

            if len(sel) != 1:
                raise ValueError(f"Found {len(sel)} atoms for residue {res} atom name {atom_name}")

            row[k] = sel.positions[0]

            if sel in masked_atoms:
                if k == 0:
                    seq[0].append(aa)
                    seq[1].append(aa)
                continue
            else:
                if k == 0:
                    seq[0].append(aa)
                    seq[1].append("<mask>")
                    
        coords.append(row)

    return np.stack(coords, axis=0), seq  


def load_esm_model(device: torch.device):
    """Load ESM-IF1 model and alphabet onto `device`."""
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval().to(device)
    return model, alphabet


@torch.no_grad()
def run_inverse_folding(
    model,
    alphabet,
    coords: np.ndarray,  
    partial_seq: list[str],
    device: torch.device,
    temperature: float = 1.0,
) -> str:
    """
    Run ESM-IF1 on one structure.

    Returns a dict with:
        sequence : sampled sequence string (length L), drawn from the
                   temperature-scaled per-position distribution
        logits   : per-position log-probabilities for the 20 standard AAs
                   shape (L, 20) serialised as a nested list
    """

    batch_converter = esm.inverse_folding.util.CoordBatchConverter(alphabet)

    L = len(coords)
    # Convert to batch format
    batch_coords, confidence, _, _, padding_mask = (
            batch_converter([(coords, None, None)], device=device)
        )

    # Start with prepend token
    mask_idx = model.decoder.dictionary.get_idx('<mask>')
    sampled_tokens = torch.full((1, 1+L), mask_idx, dtype=int)
    sampled_tokens[0, 0] = model.decoder.dictionary.get_idx('<cath>')

    if partial_seq is not None:
        for i, c in enumerate(partial_seq):
            sampled_tokens[0, i+1] = model.decoder.dictionary.get_idx(c)

    # Save incremental states for faster sampling
    incremental_state = dict()

    # Run encoder only once
    encoder_out = model.encoder(batch_coords, padding_mask, confidence)

    # Make sure all tensors are on the same device if a GPU is present
    if device:
        sampled_tokens = sampled_tokens.to(device)

    # Decode one token at a time
    for i in range(1, L+1):
        logits, _ = model.decoder(
            sampled_tokens[:, :i],
            encoder_out,
            incremental_state=incremental_state,
        )
        logits = logits[0].transpose(0, 1)

        logits /= temperature
        probs = F.softmax(logits, dim=-1)
        if sampled_tokens[0, i] == mask_idx:
            sampled_tokens[:, i] = torch.multinomial(probs, 1).squeeze(-1)

    sampled_seq = sampled_tokens[0, 1:]

    # Convert back to string via lookup
    return ''.join([model.decoder.dictionary.get_tok(a) for a in sampled_seq])


def worker(
    gpu_id: int,
    task_queue: mp.Queue,
    result_queue: mp.Queue,
    data_dir: Path,
    prefix: str,
    results_dir: Path,
    temperature: float,
    masked_atoms_selection: str,
):
    device = torch.device(f"cuda:{gpu_id}")
    print(f"[GPU {gpu_id}] Loading ESM-IF1 model …", flush=True)
    model, alphabet = load_esm_model(device)
    print(f"[GPU {gpu_id}] Model ready.", flush=True)

    while True:
        item = task_queue.get()
        if item is None:
            break

        state_i, struct_j = item
        pdb_path = data_dir / f"hidden_state_{state_i}" / f"{prefix}{struct_j}.pdb"
        out_dir = results_dir / f"state_{state_i}"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{prefix}{struct_j}.json"

        if out_path.exists():
            result_queue.put(("skip", state_i, struct_j, gpu_id))
            continue
        try:
            # Load (masked) coordinates
            coords, (seq, masked_seq) = load_backbone_coords(pdb_path, masked_atoms=masked_atoms_selection)
            # Inference!
            sequence = run_inverse_folding(model, alphabet, coords, masked_seq, device, temperature=temperature,)
            # Save
            payload = {
                "state": state_i,
                "structure": struct_j,
                "sequence": sequence,
            }
            out_path.write_text(json.dumps(payload))
            result_queue.put(("ok", state_i, struct_j, gpu_id))

        except Exception as e:
            err = traceback.format_exc()
            result_queue.put(("error", state_i, struct_j, gpu_id, str(e), err))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ESM-IF1 inverse folding on structure ensembles."
    )

    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Directory containing hidden_state_X folders",
    )

    parser.add_argument(
        "--results-dir",
        type=Path,
        default="results",
        help="Directory for output json files",
    )

    parser.add_argument("--file-prefix",
                        type=str,
                        default="structure",
                        help="The prefix of the structure file. The full file name is composed of {file_prefix}_{index}.pdb")

    parser.add_argument(
        "--num-states",
        type=int,
        default=8,
    )

    parser.add_argument(
        "--num-structures",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "--workers-per-gpu",
        type=int,
        default=1,
        help="Number of worker processes per GPU",
    )

    parser.add_argument(
        "--temperature",
        type=float,
        default=1.0,
        help="Sampling temperature for inverse folding",
    )

    parser.add_argument(
        "--masked-atoms-selection",
        type=str,
        default=None,
        help='MDAnalysis atom selection string, e.g. "protein and resid 94"',
    )

    return parser.parse_args()


def main():

    num_gpus = torch.cuda.device_count()

    if num_gpus < 1:
        raise RuntimeError("No CUDA GPUs detected; NUM_GPUS == 0")

    args = parse_args()

    args.results_dir.mkdir(parents=True, exist_ok=True)

    num_workers = num_gpus * args.workers_per_gpu

    mp.set_start_method("spawn", force=True)
    ctx = mp.get_context("spawn")

    tasks = [
        (i, j)
        for i in range(1, args.num_states + 1)
        for j in range(1, args.num_structures + 1)
    ]

    total = len(tasks)

    print(f"Detected GPUs           : {num_gpus}")
    print(f"Workers per GPU         : {args.workers_per_gpu}")
    print(f"Total workers           : {num_workers}")
    print(f"Temperature             : {args.temperature}")
    print(f"Masked atom selection   : {args.masked_atoms_selection}")
    print(f"Total tasks             : {total}  ({args.num_states} states × {args.num_structures} structures)")
    
    task_queue: mp.Queue = ctx.Queue()
    result_queue: mp.Queue = ctx.Queue()

    for task in tasks:
        task_queue.put(task)
    
    workers = []
    
    for worker_id in range(num_workers):
        gpu_id = worker_id % num_gpus
        p = ctx.Process(
            target=worker,
            args=(gpu_id,
                  task_queue, 
                  result_queue, 
                  args.data_dir,
                  args.file_prefix,
                  args.results_dir,
                  args.temperature, 
                  args.masked_atoms_selection,),
            daemon=True
                       )
        p.start()
        workers.append(p)
        task_queue.put(None)
    
    done = 0
    errors = []

    while done < total:

        item = result_queue.get()

        status, state_i, struct_j, gpu_id = (
            item[0],
            item[1],
            item[2],
            item[3],
        )

        done += 1

        if status == "ok":
            if done % 100 == 0:
                print(
                    f"[{done}/{total}] GPU {gpu_id} ✓ "
                    f"state={state_i} struct={struct_j}",
                    flush=True,
                )

        elif status == "skip":
            if done % 100 == 0:
                print(
                    f"[{done}/{total}] GPU {gpu_id} — skip "
                    f"state={state_i} struct={struct_j}",
                    flush=True,
                )

        elif status == "error":
            msg = item[4]

            print(
                f"[{done}/{total}] GPU {gpu_id} ✗ "
                f"state={state_i} struct={struct_j} ERROR: {msg}",
                flush=True,
            )

            errors.append((state_i, struct_j, item[5]))

    for p in workers:
        p.join()

    print(f"\nFinished. {total - len(errors)} succeeded, {len(errors)} failed.")


if __name__ == "__main__":
    main()
