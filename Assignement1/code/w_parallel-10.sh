#!/bin/bash
#PBS -q dssc
#PBS -l nodes=1:ppn=48
#PBS -l walltime=12:00:00


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="10"
LAUNCHES=3

for i in $(seq ${LAUNCHES})
do
	/usr/bin/time mpirun  --mca btl '^openib' -np  4 mpi_pi.x  40000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np  8 mpi_pi.x  80000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 12 mpi_pi.x 120000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 16 mpi_pi.x 160000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 20 mpi_pi.x 200000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 24 mpi_pi.x 240000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 28 mpi_pi.x 280000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 32 mpi_pi.x 320000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 36 mpi_pi.x 360000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 40 mpi_pi.x 400000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 44 mpi_pi.x 440000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 48 mpi_pi.x 480000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
done
echo "Done."
