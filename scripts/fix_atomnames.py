#!/usr/bin/env python3
import sys
import os

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


def main():
    if len(sys.argv) < 2:
        print("Usage: rename_ile_cd.py [-inplace] file1.pdb [file2.pdb ...]")
        sys.exit(1)

    inplace = False
    files = []

    for arg in sys.argv[1:]:
        if arg == "-inplace":
            inplace = True
        else:
            files.append(arg)

    for pdb in files:
        if not os.path.isfile(pdb):
            print(f"Skipping non-existent file: {pdb}", file=sys.stderr)
            continue
        process_file(pdb, inplace=inplace)


if __name__ == "__main__":
    main()
