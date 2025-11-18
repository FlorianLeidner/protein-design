import os
import argparse

def rename_cd_to_cd1(line: str) -> str:
    """Rename atom CD to CD1 for ILE residues (preserving PDB format)."""
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return line

    resname = line[17:20]
    atomname = line[12:16].strip()

    if resname == "ILE" and atomname == "CD":
        # Replace columns 13–16 (Python indices 12–16) with right-justified 'CD1'
        return f"{line[:12]}{'CD1':>4}{line[16:]}"
    else:
        return line


def process_file(filename: str, inplace: bool = False):
    """Process a PDB file and rename CD to CD1 in ILE residues."""
    with open(filename, "r") as f:
        lines = f.readlines()

    modified = [rename_cd_to_cd1(line) for line in lines]

    if inplace:
        backup = filename + ".bak"
        os.replace(filename, backup)
        with open(filename, "w") as f:
            f.writelines(modified)
        print(f"Edited in-place: {filename} (backup: {backup})")
    else:
        outname = filename.rsplit(".", 1)[0] + ".renamed.pdb"
        with open(outname, "w") as f:
            f.writelines(modified)
        print(f"Wrote: {outname}")

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(description="Make sure atom names are rfdiffusion conform")

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

    args = parser.parse_args()

    for filename in args.filenames:
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} doe not exist")

    return args


def main():

    args = parse_args()

    for pdb in args.filenames:
        process_file(pdb, inplace=args.inplace)


if __name__ == "__main__":
    main()
