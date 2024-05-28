#!/bin/bash
#SBATCH --job-name=4p-32Df
#SBATCH --partition=thcp1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=out/out_%j.out
#SBATCH --error=err/err_%j.err
export OMP_NUM_THREADS=48
yhrun  ~/fcl/proton_32D/build/qlat.x
