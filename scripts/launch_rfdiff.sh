#!/bin/bash -e

exclude=$(scontrol show node | awk '/NodeName=/{n=$1} /Gres=/{if($0 ~ /nvidia-geforce-gtx-1080/ || $0 ~ /nvidia-geforce-gtx-1080-ti/) {sub("NodeName=","",n); nodes=nodes","n}} END {gsub(/^,/,"",nodes); print nodes}')

# Add nodes that have issues with rfdiffusion/singularity
exclude="${exclude},n71-[32,34],n72-[02,12,30,40],n73-[22,30],n74-40"

for i in {1..8}; do
        INDIR="${HOME}/CaADH/analysis/data/msm/structures/2025_09_14/hidden_state_${i}"
        OUTDIR="${HOME}/CaADH/design/rfdiff/output/hidden_state_${i}"
        sbatch -J "HS${i}" --export=ALL,INDIR=$INDIR,OUTDIR=$OUTDIR --exclude=$exclude jobscript_rfdiff.sbatch
done
