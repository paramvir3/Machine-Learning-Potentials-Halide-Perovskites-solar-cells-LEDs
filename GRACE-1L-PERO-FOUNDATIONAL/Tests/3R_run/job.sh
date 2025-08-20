#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=k2-epsrc-gpu-a100
#SBATCH --time=10:00:00
#SBATCH --mem=40GB
#SBATCH --gres=gpu:a100:1

module purge
module load libs/nvidia-cuda/12.8.0/bin
module load compilers/gcc/14.1.0
module load apps/anaconda3/2024.06/bin
module load mpi/openmpi/5.0.3/gcc-14.1.0

#module load libs/blas/3.12.0/gcc-14.1.0
#module load lapack/3.12.0/gcc-9.3.0

source activate /mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/grace

source /mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/lammps-new/plumed2/sourceme.sh

lmp=/mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/grace/lammps/build_v100/lmp

#export CUDA_VISIBLE_DEVICES=$((OMPI_COMM_WORLD_RANK % 2))

#TF_CPP_MIN_LOG_LEVEL=1 mpirun -np 2 --bind-to none --oversubscribe  $lmp -in in.lammps

$lmp -in in.lammps > output 
