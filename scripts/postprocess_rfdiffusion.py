import os
import sys
import json
import pickle
import warnings
import itertools

import argparse

import numpy as np
import mdtraj as md


from scipy.spatial import distance_matrix

# Outfile header
header_template = ("index,name,inpdb,intrb,overlap,major_clashes,minor_clashes,ligand_hbonds,ligand_contacts,"
                   "ligand_sasa (nm**2),pass_filter\n")
# logic functions used with filter
logic_functions = {
    "eq": lambda x, ref: x==ref,
    "gt": lambda x, ref: x>ref,
    "lt": lambda x, ref: x<ref,
    "lt_eq": lambda x, ref: x <= ref,
    "gt_eq": lambda x, ref: x >= ref
}


def calculate_distance_matrix(coords: np.ndarray):
    dmat = distance_matrix(coords, coords)
    row_indices, column_indices = np.tril_indices(dmat.shape[0], k=-1)
    dmat_tril = dmat[row_indices, column_indices]

    return dmat_tril, row_indices, column_indices

def calculate_ligand_sasa(system, ligand) -> float:

    ligand_atom_ids = system.topology.select(ligand)
    sasa = md.shrake_rupley(system, mode="atom")[0]

    return np.sum(sasa[ligand_atom_ids])

def calculate_ligand_hbonds(system, ligand) -> int:
    topology = system.topology
    ligand_ids = topology.select(ligand)

    hbonds = md.baker_hubbard(system)

    ligand_hbonds = 0

    for hbond in hbonds: # donor index, hydrogen index and acceptor index
        is_ligand_atom = [i in ligand_ids for i in hbond]
        if any(is_ligand_atom) and not all(is_ligand_atom):
            ligand_hbonds += 1

    return ligand_hbonds

def calculate_ligand_contacts(system, ligand, contact_cutoff=0.4) -> int:
    topology = system.topology
    ligand_ids = topology.select(ligand)
    ligand_resids = np.unique([topology.atom(i).residue.index for i in ligand_ids])
    protein_resids = [res.index for res in topology.residues if res.index not in ligand_resids]

    index_pairs = list(itertools.product(ligand_resids, protein_resids))

    distances = md.compute_contacts(system, index_pairs, scheme="closest-heavy")
    n_contacts = np.sum(distances[0]<contact_cutoff)
    return n_contacts



def check_distances(system, overlap_cuttoff: float = 0.01,
                    error_cutoff: float = 0.1, warning_cutoff: float = 0.15, verbose: bool = False):
    """
    Check for atom overlaps and clashes
    :param system:
    :param overlap_cuttoff:
    :param error_cutoff:
    :param warning_cutoff:
    :param verbose:
    :return:
    """

    overlap = 0
    major_clash = 0
    minor_clash = 0

    dmat, row_indices, column_indices= calculate_distance_matrix(system.xyz[0])

    topology = system.topology

    G = topology.to_bondgraph()

    for r, i, j in zip(dmat, row_indices, column_indices):
        atm1, atm2 = topology.atom(i), topology.atom(j)
        if G.has_edge(atm1, atm2) or G.has_edge(atm2, atm1):
            continue
        if r < overlap_cuttoff:
            overlap += 1
            if verbose:
                warnings.warn(f"Atoms {atm1} and {atm2} are only {r} angstrom apart", Warning)
        elif r < error_cutoff:
            major_clash += 1
            if verbose:
                warnings.warn(f"Atoms {atm1} and {atm2} are only {r} angstrom apart", Warning)
        elif r < warning_cutoff:
            minor_clash += 1
            if verbose:
                warnings.warn(f"Atoms {atm1} and {atm2} are only {r} angstrom apart",Warning)

    return overlap, major_clash, minor_clash


def check_ligand(system, ligand: str):

    topology = system.topology
    ligand_atom_ids = topology.select(ligand)

    if not len(ligand_atom_ids):
        warnings.warn(f"'{ligand}' does not match any atoms", Warning)

def mask_designed_res(system, mask: list[bool]):

    topology = system.topology

    atom_indices = []

    i = 0
    for res in topology.residues:
        if res.is_protein:
            for atm in res.atoms:
                if atm.is_sidechain and not mask[i]:
                    continue
                else:
                    atom_indices.append(atm.index)
            i += 1
        else:
            atom_indices.extend([atm.index for atm in res.atoms])

    return system.atom_slice(atom_indices)

def validate_filter(filename):

    valid_keywords = {"overlap", "major_clash", "minor_clash",
                      "ligand_hbonds", "ligand_contacts", "ligand_sasa"}

    docstring = ("\n".join(["You can specify filters in a json file. The file format is:",
                            "[[KEY, LOGIC, VALUE], ...]",
                            "KEY specifies the observable, valid values are:\n",
                            "* "+"\n* ".join(valid_keywords),
                            "\n",
                            "LOGIC specifies a logic function:\n",
                            "* "+"\n* ".join(list(logic_functions.keys())),
                            "\n",
                            "VALUE specifies the target value\n"
                            ]))

    if not filename:
        print(docstring)
        sys.exit(0)

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"No such file {filename}")
    with open(filename, "rb") as fh:
        filters = json.load(fh)

    for (keyword, logic, value) in filters:
        if keyword not in valid_keywords:
            raise ValueError(f"{keyword} is not a valid keyword.")
        for k in logic.split("_"):
            if k.lower() not in logic_functions:
                raise ValueError("{k} not a valid logic function")
        if not any([isinstance(value, float), isinstance(value, int)]):
            raise TypeError(f"The filter value must be float or int, got {type(value)} instead")

def apply_filters(filters, data_dict):

    tests = []

    for (keyword, logic, ref_value) in filters:
        if keyword not in data_dict: # This can happen when structure is portein only but ligand keywords are provided
            warnings.warn(f"{keyword} was not calculated for {data_dict['filename']}", Warning)
            return False
        else:
            tests.append(logic_functions[logic](data_dict[keyword], ref_value))

    return all(tests)

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(description="Quality check for RFDifussion output files")

    parser.add_argument("-f",
                        "--filename",
                        dest= "filename",
                        nargs= "+",
                        metavar= "FILE",
                        help= "One or more pdb files")

    parser.add_argument("--ligand",
                        dest= "ligand",
                        nargs= "+",
                        default= None,
                        type= str,
                        metavar= "STRING",
                        help= "Ligand atom selection. If multiple input files with ligands are provided the ligand\
                         selection must match the number of input files."
                        )

    parser.add_argument("-o",
                        "--outfile",
                        dest= "outfile",
                        default= "postproc.log",
                        type= str,
                        metavar= "FILE",
                        help="Write data to outfile. Can write to an existing file with the -a flag."
                        )

    parser.add_argument("-a",
                        "--append",
                        dest= "append",
                        default= False,
                        action= "store_true",
                        help= "Provide this flag if you want to append to an existing file")

    parser.add_argument("--filter",
                        dest= "filter",
                        default= None,
                        type= str,
                        metavar= "FILE",
                        help= "Filter structures by the following parameter. Use this flag with an empty string to get all options"
                        )

    args  = parser.parse_args()

    if args.filter is not None:
        validate_filter(args.filter)

    args.n_files = len(args.filename)

    for filename in args.filename:
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} doe not exist")

    if args.ligand is not None and args.n_files > 1:
        if len(args.ligand) != args.n_files:
            if len(args.ligand) == 1:
                args.ligand = [args.ligand[0], ] * args.n_files
            else:
                parser.error(f"Expected {args.n_files} ligands selections but got {len(args.ligand)}")

    if os.path.isfile(args.outfile) and not args.append:
        parser.error("Outfile already exists. Either remove existing file or provide the -a flag to append to file.")

    elif os.path.isfile(args.outfile) and args.append:
        with open(args.outfile, "r") as fh:
            header = fh.readlines()[0]
        if header != header_template:
            parser.error(f"Can only append to files that start with the header: \n {header_template}")

    return args

def main():
    args = parse_args()
    filenames = args.filename

    # Prepare score file
    outfile = args.outfile
    if not args.append:
        i0 = 0
        with open(outfile, "w") as fh:
            fh.write(header_template)
    else: # Determine index
        with open(outfile, "r") as fh:
            lines = fh.readlines()
        i0 = len(lines)-1  # Subtract header

    # Load filters if required
    if args.filter is not None:
        with open(args.filter, "rb") as fh:
            filters = json.load(fh)

    # First perform protein specific checks:
    for i, filename in enumerate(filenames, start=1):

        filename = os.path.abspath(filename)

        data_dict = {}

        base = ".".join(os.path.basename(filename).split(".")[:-1])

        trbfile = f"{'.'.join(filename.split('.')[:-1])}.trb"

        with open(trbfile, "rb") as fh:
            trb = pickle.load(fh)

        data = [f"{i0 + i}", base, filename, trbfile]
        data_dict["filename"] = filename

        system_raw = md.load(filename)

        # Mask the ghost sidechains of designed residues
        system = mask_designed_res(system_raw, trb["mask_1d"])

        clashes = check_distances(system, overlap_cuttoff=0.01, error_cutoff=0.1, warning_cutoff=0.15)
        data.extend([f"{c}" for c in clashes])
        data_dict["overlap"] = clashes[0]
        data_dict["major_clash"] = clashes[1]
        data_dict["minor_clash"] = clashes[2]

        # For Protein/Ligand systems perform additional tests
        ligands = args.ligand

        if ligands is None: # No Ligands, Nothing to do
            data.extend(["", "", ""])
        else:
            ligand = ligands[i-1]
            check_ligand(system, ligand)

            ligand_hbonds = calculate_ligand_hbonds(system, ligand)
            data.append(f"{ligand_hbonds}")
            data_dict["ligand_hbonds"] = ligand_hbonds
            ligand_contacts = calculate_ligand_contacts(system, ligand)
            data.append(f"{ligand_contacts}")
            data_dict["ligand_contacts"] = ligand_contacts
            # Check if we have any overlaps which would cause the sasa algorithm to fail
            if clashes[0] > 0:
                data.append("")
            else:
                ligand_sasa = calculate_ligand_sasa(system, ligand)
                data.append(f"{ligand_sasa:.2f}")
                data_dict["ligand_sasa"] = ligand_sasa

        if args.filter is None:
            data.append("True")
        else:
            pass_filter = apply_filters(filters, data_dict)
            data.append(str(pass_filter))

        with open(outfile, "a") as fh:
            fh.write(", ".join(data)+"\n")

if __name__ == "__main__":
    main()