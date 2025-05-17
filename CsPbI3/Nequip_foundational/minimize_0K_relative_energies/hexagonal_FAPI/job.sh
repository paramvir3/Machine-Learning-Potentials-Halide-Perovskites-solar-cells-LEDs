#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=k2-epsrc
#SBATCH --time=05:00:00
#SBATCH --ntasks=16
#SBATCH --mem=100GB

module purge

module load apps/anaconda3/2024.06/bin
module load libs/nvidia-cuda/12.4.0/bin
module load gcc
module load libs/gcc/14.1.0
module load mpi/openmpi/5.0.3/gcc-14.1.0

export PATH="/users/pahlawat/gridware/share/python/3.10.5/bin:$PATH"
export PYTHONPATH="/users/pahlawat/gridware/share/python/3.10.5/lib/python3.9/site-packages/:$PYTHONPATH"

source activate /mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/multi-gpu-7net
#
python3 relax.py >> output
