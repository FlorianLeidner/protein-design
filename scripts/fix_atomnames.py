import os
import argparse

import MDAnalysis as mda

def process_file(filename: str, inplace: bool = False, verbose: bool = False):

    u = mda.Universe(filename)
    protein = u.select_atoms("protein")


    # Fix atom names
    for atm in protein:
        if atm.name == "OT1":
            atm.name = "O"
            if verbose:
                print(f"Renamed terminal oxygen 'OT1' to 'O' Residue: {atm.resid}:{atm.resname}")
        elif atm.name == "OT2":
            atm.name = "OXT"
            if verbose:
                print(f"Renamed terminal oxygen 'OT2' to 'OXT' Residue: {atm.resid}:{atm.resname}")

        if atm.resname == "ILE" and atm.name == "CD1":
            atm.name = "CD"
            if verbose:
                print(f"Renamed ILE delta carbon from 'CD1' to 'CD' Residue: {atm.resid}:{atm.resname}")

    # Fix record type
    hetatms = u.select_atoms("all and not (protein or nucleic)")
    for atm in hetatms:
        atm.record_type = "HETATM"

    if inplace:
        u.atoms.write(filename)
        if verbose:
            print(f"Edited in-place: {filename}")
    else:
        outname = filename.rsplit(".", 1)[0] + ".renamed.pdb"
        u.atoms.write(outname)
        if verbose:
            print(f"Wrote: {outname}")

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(description="Neural Networks are picky eaters. This script fixes atom names that "
                                                 "can cause the catastrophic failure of RFDiffusion and HPacker")

    parser.add_argument("-f",
                        "--filenames",
                        dest= "filenames",
                        nargs= "+",
                        metavar= "FILE",
                        help= "One or more pdb files")

    parser.add_argument("--inplace",
                        dest="inplace",
                        action="store_true",
                        default=False,
                        help="Modify files in place. If False (default), the input file will be backed up.")

    parser.add_argument("-v",
                        dest = "verbose",
                        action = "store_true",
                        default= False,
                        help= "Make some noise")

    args = parser.parse_args()

    for filename in args.filenames:
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} doe not exist")

    return args


def main():

    args = parse_args()

    for pdb in args.filenames:
        process_file(pdb, inplace=args.inplace, verbose=args.verbose)


if __name__ == "__main__":
    main()
