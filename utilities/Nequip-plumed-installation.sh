#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=k2-epsrc-gpu
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=100GB
#SBATCH --gres=gpu:a100:1

module purge
module load libs/gcc/14.1.0
module load libs/nvidia-cuda/11.8.0/bin
module load mpi/openmpi/5.0.3/gcc-14.1.0
module load libs/blas/3.12.0/gcc-14.1.0
module load lapack/3.12.0/gcc-9.3.0
module load apps/anaconda3/2024.06/bin

source activate /mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/lammps-new/lammps_latest_nequip/neqlmp
source /mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/lammps-new/plumed2/sourceme.sh
plumed_dir="/mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/lammps-new/plumed2/"

git clone --depth=1 https://github.com/lammps/lammps
git clone https://github.com/mir-group/pair_nequip
cd pair_nequip
./patch_lammps.sh ../lammps/
cd ..

cd lammps
rm -rf build_a100
mkdir build_a100
cd build_a100

cmake \
    -D CMAKE_BUILD_TYPE=Release \
    -D LAMMPS_EXCEPTIONS=yes \
    -D PKG_MANYBODY=yes \
    -D PKG_EXTRA-FIX=yes \
    -D PKG_EXTRA-PAIR=yes \
    -D PKG_EXTRA-DUMP=yes \
    -D PKG_MOLECULE=yes \
    -D Python_EXECUTABLE=$(which python3) \
    -D CMAKE_PREFIX_PATH="`python -c 'import torch;print(torch.utils.cmake_prefix_path)'`;${plumed_dir}" \
    -D PKG_PLUMED=yes  -D PLUMED_MODE=shared  -D DOWNLOAD_PLUMED=no \
    -D MKL_INCLUDE_DIR=`python -c "import sysconfig;from pathlib import Path;print(Path(sysconfig.get_paths()[\"include\"]).parent)"` \
    -D CUDA_TOOLKIT_ROOT_DIR=/opt/gridware/depots/54e7fb3c/el8/pkg/libs/nvidia-cuda/11.8.0/bin \
    -D CUDNN_LIBRARY=/mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/lammps-new/lammps_latest_nequip/cudnn/lib  \
    -D CUDNN_INCLUDE_DIR=/mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/lammps-new/lammps_latest_nequip/cudnn/include \
    ../cmake
make -j4
