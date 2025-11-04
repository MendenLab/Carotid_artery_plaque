#!/bin/bash
#SBATCH -o run_velocyto.log
#SBATCH -e run_velocyto.log
#SBATCH --job-name=JCPU_s8_4
#SBATCH --ntasks=1
#SBATCH -c 16
#SBATCH --partition=cpu_p
#SBATCH --mem-per-cpu=16G
#SBATCH --qos=cpu_normal
##SBATCH --gres=gpu:2
#SBATCH --time=14:00:00
#SBATCH --nice=10000

set -eo pipefail

source $HOME/.bashrc

echo "Starting jupyter..."

#chmod 600 $HOME/slurm_jupyter_$SLURM_JOB_ID.job
#chmod 600 ./slurm_jupyter_$SLURM_JOB_ID.job

#eval "$(conda shell.bash hook)"
conda deactivate
conda activate velocyto

sh run_velocyto.sh
#conda activate scarches
#jupyter-lab --port 8888 --no-browser --ip=0.0.0.0 --notebook-dir=$HOME/
#jupyter-lab --notebook-dir=$HOME/ --no-browser  --ip=0.0.0.0
