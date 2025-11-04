#!/bin/bash
#SBATCH --job-name=pyscenic
#SBATCH -o pyscenic_output/grn_%x_%j.out
#SBATCH -e pyscenic_output/grn_%x_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --partition=sapphire
#SBATCH --time=1-00:00:00
#SBATCH --nice=10000

set -euo pipefail



echo "Loading anaconda"
module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"

echo "Activate conda env"
conda activate pyscenic

out_fn='../../data/GRCh38-p14-Gencode_v44_with_intron_adj2.csv'
f_loom_path_scenic='../../data/GRCh38-p14-Gencode_v44_with_intron_pyscenic.loom'

PYTHONNOUSERSITE=1 conda run --no-capture-output -n pyscenic pyscenic grn \
${f_loom_path_scenic} \
"../../pyscenic_data/allTFs_hg38.txt" \
-o ${out_fn} \
--gene_attribute Gene \
--num_workers 32