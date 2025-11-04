#!/bin/bash
#SBATCH --job-name=pyscenic
#SBATCH -o pyscenic_output/ctx_%x_%j.out
#SBATCH -e pyscenic_output/ctx_%x_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --partition=bigmem
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

echo $(which pyscenic)

echo "Python in env:"
which python
python -V

adj_fn='../../data/GRCh38-p14-Gencode_v44_with_intron_adj2.csv'
f_loom_path_scenic='../../data/GRCh38-p14-Gencode_v44_with_intron_pyscenic.loom'
reg_fn='../../data/GRCh38-p14-Gencode_v44_with_intron_reg.csv'


# path to ranking databases in feather format
f_db_glob="../../pyscenic_data/v10/*feather"
f_db_names=$(echo $f_db_glob)   # literal glob string
f_db_names=$(ls ../../pyscenic_data/v10/*feather 2>/dev/null | tr '\n' ' ')
#echo ${f_db_names}


# path to motif databases
#f_motif_path="../../pyscenic_data/v9/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
f_motif_path="../../pyscenic_data/v10/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

##PYTHONNOUSERSITE=1 conda run --no-capture-output -n pyscenic pyscenic grn ${f_loom_path_scenic} "../../data/allTFs_hg38.txt" -o ${out_fn} --#num_workers 44


echo "Running pyscenic ctx..."
#PYTHONNOUSERSITE=1 conda run --no-capture-output -n pyscenic pyscenic ctx ${adj_fn} ${f_db_names} --annotations_fname ${f_motif_path} --expression_mtx_fname #${f_loom_path_scenic} --output ${reg_fn} --all_modules --num_workers 32 # --no_pruning --transpose 

PYTHONNOUSERSITE=1 conda run --no-capture-output -n pyscenic pyscenic ctx \
  ${adj_fn} ${f_db_names} \
  --annotations_fname ${f_motif_path} \
  --expression_mtx_fname ${f_loom_path_scenic} \
  --output ${reg_fn} \
  --all_modules \
  --rank_threshold 10000 \
  --auc_threshold 0.05 \
  --nes_threshold 2.5 \
  --max_similarity_fdr 0.01 \
  --top_n_targets 300 \
  --min_genes 10 \
  --num_workers 32
