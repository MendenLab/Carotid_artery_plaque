#!/bin/bash
#SBATCH --job-name=pyscenic
#SBATCH -o pyscenic_output/aucell_%x_%j.out
#SBATCH -e pyscenic_output/aucell_%x_%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --partition=sapphire
#SBATCH --time=1-00:00:00
#SBATCH --nice=10000

set -euo pipefail


echo "Loading anaconda"
module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"

echo "Activate conda env"
PYS_ENV="/data/gpfs/projects/punim2121/anaconda3/envs/pyscenic"
#conda activate "$PYS_ENV"
conda activate pyscenic


f_loom_path_scenic='../../data/GRCh38-p14-Gencode_v44_with_intron_pyscenic.loom'
reg_fn='../../data/GRCh38-p14-Gencode_v44_with_intron_reg.csv'
aucell_output_fn='../../data/GRCh38-p14-Gencode_v44_with_intron_aucell_output_fn.loom'



#pyscenic aucell ${f_loom_path_scenic} ${reg_fn} --output ${aucell_output_fn} --num_workers 4


echo "Running pyscenic aucell..."
PYTHONNOUSERSITE=1 conda run --no-capture-output -n pyscenic pyscenic aucell ${f_loom_path_scenic} ${reg_fn} --output ${aucell_output_fn} --num_workers 4