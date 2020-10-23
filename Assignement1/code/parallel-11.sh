#!/bin/bash
#PBS -q dssc
#PBS -l nodes=1:ppn=48
#PBS -l walltime=06:15:00


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="100000000000"
LAUNCHES=3

for procs in 1 4 8 12 16 20 24 28 32 36 40 44 48 ; do
	echo "Parallel execution. P= ", ${procs}
	echo "N= " ${MOVES}
	for i in $(seq ${LAUNCHES})
	do
		/usr/bin/time mpirun  --mca btl '^openib' -np ${procs} mpi_pi.x  ${MOVES} >>parallel-${MOVES}.output 2>>parallel-${MOVES}.error
	done
done
echo "Done."
