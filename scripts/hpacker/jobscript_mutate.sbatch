#!/bin/bash -e
#SBATCH -p p16,p20,p32,p40
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

source ${HOME}/software/miniconda/etc/profile.d/conda.sh
conda activate hpacker


script="${HOME}/repos/protein-design/scripts/hpacker/mutate.py"

if [ -z "$MUTATION" ]; then
    echo "Provide a valid instruction string for HPacker with the export=ALL,MUTATION=XXX statement"
    exit 1
fi

echo "${MUTATION}"

python $script "${MUTATION}"
 
