#!/bin/bash
#PBS -q dssc
#PBS -l nodes=1:ppn=48
#PBS -l walltime=00:10:00


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="10"

/usr/bin/time mpirun  --mca btl '^openib' -np 48 mpi_pi.x 480000000000 >>r_w_parallel-${MOVES}.output 2>>r_w_parallel-${MOVES}.error
echo "Done."
