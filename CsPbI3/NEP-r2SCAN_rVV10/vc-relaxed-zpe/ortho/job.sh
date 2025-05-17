#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=k2-epsrc-gpu-v100
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --gres=gpu:v100:1

module purge
module load libs/gcc/14.1.0
module load libs/nvidia-cuda/12.4.0/bin

export PATH="/opt/gridware/depots/54e7fb3c/el8/pkg/libs/nvidia-cuda/12.4.0/bin/bin/:$PATH"
export LD_LIBRARY_PATH="/opt/gridware/depots/54e7fb3c/el8/pkg/libs/nvidia-cuda/12.4.0/bin/lib64/:$LD_LIBRARY_PATH"

/mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/GPUMD_new_a100/GPUMD/src/gpumd < run.in &> output

