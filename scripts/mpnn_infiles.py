import os
import json
import pickle
import warnings
import argparse

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def make_resmask(trb: dict) -> str:

    reslist = []
    for i, (chain, resnum) in enumerate(trb["con_hal_pdb_idx"]):
        if trb["mask_1d"]:
            reslist.append(f"{chain}{resnum}")

    return " ".join(reslist)

def bias_from_seq(seq: SeqRecord, trb: dict, bias_factor: float) -> dict:
    """
    Generate bias_dict based on sequence
    :param seq:
    :param trb:
    :param bias_factor: How much to bias towards a specific AA.
     This affects the amino acid likelihood in the following way: softmax((logit+bias)/T)
     where logit is the output of the last linear layer, T is the temperature, and bias is the aa bias_factor.
    :return:
    """

    bias_dict = {}

    chain_id = trb["con_hal_pdb_idx"][0][0]

    if not all([c == chain_id for (c, r) in trb["con_hal_pdb_idx"]]):
        raise NotImplementedError("Currently only single chain designs are supported.")

    n_res = len(trb["mask_1d"])

    if len(seq.seq) != n_res:
        raise ValueError(f"Sequence was {len(seq.seq)} but designed protein was {n_res} residues.")

    for i in range(n_res):
        if not trb["mask_1d"][i]:
            aa = seq.seq[i]
            resid = f"{chain_id}{i + 1}"

            if resid in bias_dict:
                bias_dict[resid][aa] = bias_factor
            else:
                bias_dict[resid] = {aa: bias_factor}

    return bias_dict

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(description="Prepare infiles for MPNN")

    parser.add_argument("-i",
                        dest="infile",
                        metavar="FILE",
                        required=True,
                        help="A csv file created by postprocess_rfdiffusion.py ")

    parser.add_argument("--prefix",
                        dest="prefix",
                        type=str,
                        help="prefix for mpnn files")

    parser.add_argument("--seq",
                        dest="seq",
                        default=None,
                        type=str,
                        metavar="FILE",
                        help="A sequence in fasta format required for generating bias from sequence")

    parser.add_argument("--seq_bias",
                        dest="seq_bias",
                        default=0.,
                        type=float,
                        metavar="FLOAT",
                        help="How much to bias towards a specific sequence")

    args = parser.parse_args()

    if not os.path.isfile(args.infile):
        raise FileNotFoundError(f"{args.infile} doe not exist")

    if args.seq is not None:
        if os.path.isfile(args.seq):
            try:
                args.seq = SeqIO.read(args.seq, format="fasta")
            except Exception as e:
                print(e)
                raise IOError("If target is file, it must in in fasta format")
        else:
            seq = args.seq.upper()
            seq = seq.strip()
            if not all([symbol.isalpha() for symbol in seq]):
                raise IOError(f"Sequence contains unknown characters: {seq}")
            else:
                args.seq = SeqRecord(Seq(seq), id="query")


    return args

def main():

    mpnn_mask = {}
    mpnn_input = {}
    mpnn_bias = {}

    args = parse_args()

    data = pd.read_csv(args.infile)
    data = data.astype({"name": str, "inpdb": str, "intrb": str, "major_clashes": int, "minor_clashes": int,
                        "ligand_hbonds": int,  "ligand_contacts": int, "ligand_sasa (nm**2)": float,
                        "pass_filter": bool})

    data = data[data["pass_filter"]]

    for index in data.index:

        pdbfile = data.loc[index, "inpdb"]
        mpnn_input[pdbfile] = ""

        # Mask fixed residues
        trbfile = data.loc[index, "intrb"]

        with open(trbfile, "rb") as fh:
            trb = pickle.load(fh)

        mask = make_resmask(trb)

        if mask:
            mpnn_mask[pdbfile] = mask

        if args.seq is not None:
            bias_dict = bias_from_seq(args.seq, trb)
            mpnn_bias[pdbfile] = bias_dict


    # --- Check if input files already exist --- #

    mpnn_input_file = f"{args.prefix}_input.json"
    if os.path.isfile(mpnn_input_file):
        warnings.warn(f"{mpnn_input_file} already exists. Created backup of old file", Warning)
        os.rename(mpnn_input_file, f"{mpnn_input_file}.bkp")

    with open(mpnn_input_file, "w") as fh:
        fh.write(json.dumps(mpnn_input))

    if mpnn_mask:
        mpnn_mask_file = f"{args.prefix}_mask.json"
        if os.path.isfile(mpnn_mask_file):
            warnings.warn(f"{mpnn_mask_file} already exists. Created backup of old file", Warning)
            os.rename(mpnn_mask_file, f"{mpnn_mask_file}.bkp")

        with open(mpnn_mask_file, "w") as fh:
            fh.write(json.dumps(mpnn_mask))

    if mpnn_bias:
        mpnn_bias_file = f"{args.prefix}_bias.json"
        if os.path.isfile(mpnn_bias_file):
            warnings.warn(f"{mpnn_bias_file} already exists. Created backup of old file", Warning)
            os.rename(mpnn_bias_file, f"{mpnn_bias_file}.bkp")

        with open(mpnn_bias_file, "w") as fh:
            fh.write(json.dumps(mpnn_bias))

if __name__=="__main__":
    main()