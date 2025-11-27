import os
import argparse
import warnings

from hpacker import HPacker

def files_from_dir(path: str =None) -> list[str]:
    if path is None:
        path = os.getcwd()
    else:
        path = os.path.abspath(path)

    files = []

    for name in os.listdir(path=path):
        isdir = os.path.isdir(os.path.join(path, name))
        if isdir:
            files.extend(files_from_dir(path=os.path.join(path, name)))
        else:
            suffix = name.split(".")[-1]
            if suffix == "pdb":
                files.append(os.path.join(path, name))
            else:
                continue
    return files

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(description="Mutate and repack a set of protein structures with HPacker")

    parser.add_argument("-i",
                        dest= "infile",
                        nargs= "+",
                        default=None,
                        metavar= "FILE",
                        help= "Process all listed pdb files.")

    parser.add_argument("-d",
                        dest="indir",
                        nargs="+",
                        default=None,
                        metavar="DIR",
                        help="Process all pdb files in the listed directories.")

    parser.add_argument("-o",
                        "--outdir",
                        dest="outdir",
                        default="./hpacker",
                        metavar="DIR",
                        help="output structures will be saved here.")

    parser.add_argument("--mutations",
                        dest="mutations",
                        nargs="+",
                        required=True,
                        help=("A list of residues to mutate. Each mutation is specified by a string of:"
                              "RESNUM:CHAIN:RESNAME . For homo-multimers each monomer has to be specified."))

    parser.add_argument("-v",
                        dest = "verbose",
                        action = "store_true",
                        default= False,
                        help= "Make some noise")

    args = parser.parse_args()

    for mut in args.mutations:
        if not len(mut.split(":")) == 3:
            parser.error(f"Invalid mutation: {mut}.Each mutation needs to be composed of: RESNUM:CHAIN:RESNAME.")

    args.pdbs = []

    if args.infile is not None:
        for filename in args.infiles:
            if not os.path.isfile(filename):
                raise FileNotFoundError(f"{filename} doe not exist")
            else:
                args.pdbs.append(filename)

    if args.indir is not None:
        for d in args.indir:
            if not os.path.isdir(d):
                raise NotADirectoryError(f"{d} is not a directory")
            else:
                args.pdbs.extend(files_from_dir(os.path.abspath(d)))

    if len(args.pdbs) == 0:
        parser.error("Provide either a list of pdb files or a list of directories containing pdb files")

    return args

def mutate(pdbfile: str, mutations: list[str], outdir="./hpacker", verbose: bool = False,
           num_refinement_iterations: int = 30, proximity_cutoff_for_refinement: float = 10.):

    filename = os.path.split(pdbfile)[-1]

    mutations_dict = {}

    for mutation in mutations:
        resnum, chain, resname = mutation.split(":")
        resnum = int(resnum)
        mutations_dict[(chain, resnum, " ")] = resname

    hpacker = HPacker(pdbfile, verbose=verbose)

    hpacker.reconstruct_sidechains(num_refinement_iterations = num_refinement_iterations,
                               proximity_cutoff_for_refinement = proximity_cutoff_for_refinement,
                               res_id_to_resname=mutations_dict)

    hpacker.write_pdb(os.path.join(outdir, filename))

def main():
    args = parse_args()

    if not os.path.isdir(args.outdir):
        warnings.warn(f"{args.outdir} is not a directory. Creating new directory.")
        # Make directory if it does not exist already
        os.makedirs(args.outdir, exist_ok=True)

    for pdbfile in args.pdbs:
        mutate(pdbfile, args.mutations, outdir=args.outdir, verbose=args.verbose)

if __name__=="__main__":
    main()