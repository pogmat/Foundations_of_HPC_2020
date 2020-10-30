#!/bin/bash
#PBS -q dssc
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:12:00


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="10000000000"
LAUNCHES=3

echo "Serial execution."
echo "N= " ${MOVES}
for i in $(seq ${LAUNCHES})
do
	/usr/bin/time ./pi.x ${MOVES} >>serial-${MOVES}.output 2>>serial-${MOVES}.error
done
echo "Done."
