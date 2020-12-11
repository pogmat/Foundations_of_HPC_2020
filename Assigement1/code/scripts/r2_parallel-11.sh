#!/bin/bash
#PBS -q dssc
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:15:00


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0


/usr/bin/time mpirun  --mca btl '^openib' -np 8 mpi_pi.x  100000000000 >>r2_parallel-11.output 2>>r2_parallel-11.error
echo "Done."
