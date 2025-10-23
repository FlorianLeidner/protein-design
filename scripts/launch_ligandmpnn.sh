#!/bin/bash -e

#Exclude all gtx-1080 nodes
exclude=$(scontrol show node | awk '/NodeName=/{n=$1} /Gres=/{if($0 ~ /nvidia-geforce-gtx-1080/ || $0 ~ /nvidia-geforce-gtx-1080-ti/) {sub("NodeName=","",n); nodes=nodes","n}} END {gsub(/^,/,"",nodes); print nodes}')

# Add nodes that have issues with rfdiffusion/singularity
exclude="${exclude},n71-[32,34],n72-[02,12,30,40],n73-[22,30],n74-40"

for i in {1..8}; do
  INDIR="${HOME}/CaADH/design/rfdiff/output/hidden_state_${i}"
  OUTDIR="${HOME}/CaADH/design/ligandmpnn/output"
  sbatch -J mpnn_hs${i} --export=ALL,INDIR=$INDIR,OUTDIR=$OUTDIR --array=1-50%25 jobscript_ligandmpnn.sh
done