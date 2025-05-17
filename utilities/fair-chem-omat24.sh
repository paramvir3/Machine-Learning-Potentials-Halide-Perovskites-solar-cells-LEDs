#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=k2-epsrc-gpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=200GB
#SBATCH --gres=gpu:v100:1

module purge
module load apps/anaconda3/2024.06/bin
module load libs/nvidia-cuda/12.4.0/bin
module load libs/gcc/14.1.0

wget https://raw.githubusercontent.com/FAIR-Chem/fairchem/main/packages/env.gpu.yml
mv env.gpu.yml environment.yml
conda env create --prefix=/mnt/scratch2/q13camb_scratch/perovskites/Balaz_project/fair-chem # own prefix

#conda init
#conda env create -f env.gpu.yml 
#echo yes|conda env create -f env.gpu.yml 

source activate /mnt/scratch2/q13camb_scratch/perovskites/Balaz_project/fair-chem

git clone -b omat24 https://github.com/FAIR-Chem/fairchem.git
cd fairchem
pip install -e packages/fairchem-core

pip uninstall pyg-lib

#pytest tests/core
