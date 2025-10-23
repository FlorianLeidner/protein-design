#!/bin/bash -e
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p p40,p32,p20,p16,p08
#SBATCH --gpus-per-task:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH --array=250-500%25

CONTAINER="${HOME}/containers/rf_se3_diffusion.sif"
RFDIFF_HOME="${HOME}/repos/rf_diffusion_all_atom"


INDEX=$(( SLURM_ARRAY_TASK_ID+1 ))

INFILE="${INDIR}/structure${INDEX}.pdb"

DESIGNS=10
START=$(( DESIGNS*SLURM_ARRAY_TASK_ID ))


singularity run --nv $CONTAINER -u ${RFDIFF_HOME}/run_inference.py \
        inference.output_prefix=${OUTDIR}/out \
        inference.input_pdb=$INFILE \
        inference.ligand="NAP" \
        inference.num_designs=$DESIGNS \
        inference.design_startnum=$START \
        inference.ckpt_path=${RFDIFF_HOME}/RFDiffusionAA_paper_weights.pt \
        contigmap.contigs=[\'A1-186,26-26,A213-251\'] \
        potentials.guiding_potentials=[\'type:ligand_ncontacts,weight:2\'] \
        potentials.guide_scale=2 \
        potentials.guide_decay="cubic" \
        diffuser.T=100 # Reduced diffusion cycles