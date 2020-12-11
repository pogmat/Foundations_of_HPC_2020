#!/bin/bash
#PBS -q dssc
#PBS -l nodes=1:ppn=48
#PBS -l walltime=01:56:00


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="09"
LAUNCHES=3

for i in $(seq ${LAUNCHES})
do
	/usr/bin/time mpirun  --mca btl '^openib' -np  4 mpi_pi.x  4000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np  8 mpi_pi.x  8000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 12 mpi_pi.x 12000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 16 mpi_pi.x 16000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 20 mpi_pi.x 20000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 24 mpi_pi.x 24000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 28 mpi_pi.x 28000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 32 mpi_pi.x 32000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 36 mpi_pi.x 36000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 40 mpi_pi.x 40000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 44 mpi_pi.x 44000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
	/usr/bin/time mpirun  --mca btl '^openib' -np 48 mpi_pi.x 48000000000 >>w_parallel-${MOVES}.output 2>>w_parallel-${MOVES}.error
done
echo "Done."
