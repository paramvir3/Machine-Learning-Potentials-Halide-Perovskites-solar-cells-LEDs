#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=k2-epsrc-gpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=300GB
#SBATCH --gres=gpu:a100:1

module purge

module load apps/anaconda3/2024.06/bin
module load libs/nvidia-cuda/12.4.0/bin
module load gcc
module load libs/gcc/14.1.0
module load  mpi/openmpi/5.0.3/gcc-14.1.0

echo yes|conda create --prefix=/mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/multi-gpu-7net

source activate /mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/multi-gpu-7net

pip install torch==2.2.2 --index-url https://download.pytorch.org/whl/cu121
#python test.py &> output  # test if pytorch has cuda 
       '''import torch  
          device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
          print(f"Using device: {device}")'''

git clone https://github.com/MDIL-SNU/SevenNet.git
cd SevenNet/
pip install .

cd .. 

git clone https://github.com/lammps/lammps.git lammps_sevenn --branch stable_2Aug2023_update3 --depth=1
sevenn_patch_lammps ./lammps_sevenn
pip install  mkl-include

cd ./lammps_sevenn
rm -rf build
mkdir build
cd build

cmake     -D CMAKE_BUILD_TYPE=Release     -D CMAKE_INSTALL_PREFIX=$(pwd)  -D BUILD_OMP=yes     -D BUILD_SHARED_LIBS=yes     -D LAMMPS_EXCEPTIONS=yes    \
          -D PKG_EXTRA-PAIR=yes     -D PKG_EXTRA-DUMP=yes     -D PKG_MOLECULE=yes     -D PKG_EXTRA-FIX=yes     -D CMAKE_PREFIX_PATH="`python3 -c 'import torch;print(torch.utils.cmake_prefix_path)'`"   \
          -D MKL_INCLUDE_DIR=`python3 -c "import sysconfig;from pathlib import Path;print(Path(sysconfig.get_paths()[\"include\"]).parent)"` -D CUDA_TOOLKIT_ROOT_DIR=/opt/gridware/depots/54e7fb3c/el8/pkg/libs/nvidia-cuda/12.4.0/bin/  \
          -D Python_EXECUTABLE=$(which python3)  ../cmake

make -j8
