#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=k2-epsrc-gpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=100GB
#SBATCH --gres=gpu:a100:1

module purge
module load apps/anaconda3/2024.06/bin
module load libs/nvidia-cuda/12.4.0/bin
module load libs/gcc/14.1.0
module load mpi/openmpi/5.0.3/gcc-14.1.0
module load libs/blas/3.12.0/gcc-14.1.0
module load lapack/3.12.0/gcc-9.3.0


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

export PATH="/users/pahlawat/gridware/share/python/3.10.5/bin:$PATH"
export PYTHONPATH="/users/pahlawat/gridware/share/python/3.10.5/lib/python3.9/site-packages/:$PYTHONPATH"

git clone https://github.com/plumed/plumed2.git -b v2.8.3 --depth 1
cd plumed2/
plumed_dir=${PWD}
./configure --enable-modules=all --prefix=${PWD}
make -j4
#make install
source ${PWD}/sourceme.sh
#echo "source ${plumed_dir}/sourceme.sh" >> ~/.bashrc
#echo `export PYTHONPATH="/home/mmm1188/Nequip_lammps/plumed2/lib/plumed/python/:$PYTHONPATH"` >> ~/.bashrc
#plumed_dir="/mnt/scratch2/q13camb_scratch/perovskites/installation_stuff/multi-gpu-7net_nvhpc_plumed/plumed2"

cd ..

git clone https://github.com/lammps/lammps.git lammps_sevenn --branch stable_2Aug2023_update3 --depth=1
sevenn_patch_lammps ./lammps_sevenn
pip install  mkl-include

cd lammps_sevenn
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
    -D CMAKE_PREFIX_PATH="`python3 -c 'import torch;print(torch.utils.cmake_prefix_path)'`;${plumed_dir}" \
    -D PKG_PLUMED=yes  -D PLUMED_MODE=shared  -D DOWNLOAD_PLUMED=no \
    -D MKL_INCLUDE_DIR=`python3 -c "import sysconfig;from pathlib import Path;print(Path(sysconfig.get_paths()[\"include\"]).parent)"` -D CUDA_TOOLKIT_ROOT_DIR=/opt/gridware/depots/54e7fb3c/el8/pkg/libs/nvidia-cuda/12.4.0/bin ../cmake

make -j4


