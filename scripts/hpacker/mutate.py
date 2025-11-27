import os
import sys
from glob import glob

from hpacker import HPacker


infiles = glob("./in/structure*.pdb")

for inpdb in infiles:

    filename = os.path.split(inpdb)[-1]

    outdir = "./hpacker"
    #Make directory if it does not exist already
    os.makedirs(outdir, exist_ok=True)


    mutations_string = sys.argv[1]
    mutations_string = mutations_string.strip("\"") 
    mutations = {}
    for mutation in mutations_string.split(","):
        resnum, chain, resname = mutation.split(":")
        resnum = int(resnum)
        mutations[(chain, resnum, " ")] = resname

    print(mutations)

    hpacker = HPacker(inpdb, verbose=True)

    hpacker.reconstruct_sidechains(num_refinement_iterations = 30, 
                               proximity_cutoff_for_refinement = 10.0, 
                               res_id_to_resname=mutations)

    hpacker.write_pdb(os.path.join(outdir, filename))
