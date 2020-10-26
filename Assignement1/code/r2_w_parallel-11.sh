#!/bin/bash
#PBS -q dssc
#PBS -l nodes=1:ppn=24
#PBS -l walltime=03:00:00


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="11"
LAUNCHES=3

for i in $(seq ${LAUNCHES})
do
	/usr/bin/time mpirun  --mca btl '^openib' -np 24 mpi_pi.x 2400000000000 >>r2_w_parallel-${MOVES}.output 2>>r2_w_parallel-${MOVES}.error
done
echo "Done."
