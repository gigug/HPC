#!/bin/bash
#PBS -l nodes=1:ppn=48
#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

mkdir MPI_strong_11

for procs in {1..48} ; do
    echo "executing on ", ${procs}, "  processors" 
    /usr/bin/time  mpirun  --mca btl '^openib' -np ${procs} ./MPI_strong 11 11 >>MPI_strong_11/MPI_strong_11_${procs}.out 2>&1
done


