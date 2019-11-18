#!/bin/bash

#SBATCH -J heat2D_v6
#SBATCH -n 16
#SBATCH -t 00:01:00
#SBATCH -o heat2D-%j.out
#SBATCH -e heat2D-%j.err
#SBATCH -D .

module purge
module load gnu/7.2.0
module load szip/gnu/2.1.1           
module load openmpi/gnu/3.0.1        
module load hdf5/gnu/openmpi/1.10.1  

mpirun -np 16 ./heat_mpi iac.h5


