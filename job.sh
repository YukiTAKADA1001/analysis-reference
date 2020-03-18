#!/bin/sh
#PBS -l select=4:ncpus=40:mpiprocs=8:ompthreads=1:jobtype=small
#PBS -l walltime=178:00:00

cd $PBS_O_WORKDIR

mpirun -np 160 pmemd.MPI -ng 8 -groupfile groupfile -rem 1
