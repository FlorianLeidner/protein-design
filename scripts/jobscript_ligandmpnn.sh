#!/bin/bash -e
#
#SBATCH -p p16,p20,p32,p40
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --time=24:00:00


SCRIPT="${HOME}/repos/protein-design/scripts/postprocess_rfdiffusion.py"

MPNN_HOME="${HOME}/repos/LigandMPNN"
CONTAINER="${HOME}/containers/ligandmpnn.sif"


DESIGNS=100

START=$(( DESIGNS*SLURM_ARRAY_TASK_ID ))
STOP=$(( START+DESIGNS))

mkdir -p "${OUTDIR}/${SLURM_ARRAY_TASK_ID}"
cd "${OUTDIR}/${SLURM_ARRAY_TASK_ID}"

# ---- Compile file list ----

PREFIX="out_"                             # File prefix
SUFFIX=".pdb"                             # File suffix

# ---- Collect existing files ----
file_list=()

for (( i=$START; i<=$STOP; i++ )); do
    filename="${INDIR}/${PREFIX}${i}${SUFFIX}"
    if [[ -f "$filename" ]]; then
        file_list+=("$filename")
    fi
done

# ---- Check if we have any files ----
if [[ ${#file_list[@]} -eq 0 ]]; then
    echo "No files found in range $START to $STOP."
    exit 1
fi


# Activate conda env
source /home/fleidne/software/miniconda/etc/profile.d/conda.sh
conda activate mdanalysis

# ---- Call the Python script ----
echo "Processing ${#file_list[@]} files..."

cat <<EOF >> filter.json
[["overlap", "eq", 0],
["major_clash", "eq", 0],
["minor_clash", "lt_eq", 10],
["ligand_sasa", "lt", 2.5],
["ligand_contacts", "gt", 20]
]
EOF

python $SCRIPT -f "${file_list[@]}" --ligand "resname NAP" -o  "${SLURM_ARRAY_TASK_ID}_rfdiff_out.csv" --prepare_mpnn_input --mpnn_prefix "${SLURM_ARRAY_TASK_ID}" --filter filter.json


# Multiple temperatures and 5 sequences per job following the recommendations in:
# https://github.com/ikalvet/heme_binder_diffusion
temperatures=(0.1 0.2 0.3)

for temp in "${temperatures[@]}"; do
    singularity run --nv $CONTAINER --model_type "ligand_mpnn" --pdb_path_multi "${SLURM_ARRAY_TASK_ID}_input.json" --fixed_residues_multi "${SLURM_ARRAY_TASK_ID}_mask.json"\
     --out_folder "${OUTDIR}/${SLURM_ARRAY_TASK_ID}/T${temp}" --checkpoint_ligand_mpnn "${MPNN_HOME}/model_params/ligandmpnn_v_32_030_25.pt" \
     --checkpoint_path_sc "${MPNN_HOME}/model_params/ligandmpnn_sc_v_32_002_16.pt" \
     --temperature $temp --number_of_batches 1 --batch_size 5 --pack_side_chains 1 --number_of_packs_per_design 1 \
      --pack_with_ligand_context 1
done
