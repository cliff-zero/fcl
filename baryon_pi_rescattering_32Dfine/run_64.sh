#!/bin/bash
#SBATCH --job-name=重楼
#SBATCH --partition=thcp1
#SBATCH --nodes=64
#SBATCH --ntasks=64
#SBATCH --output=out/out_%j.out
#SBATCH --error=err/err_%j.err
export OMP_NUM_THREADS=64
yhrun  ~/Yusheng_Gao/GYS_Contraction/baryon_pi_rescattering_32Dfine/build/qlat.x
