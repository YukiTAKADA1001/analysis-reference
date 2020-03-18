#!/bin/sh
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1:jobtype=gpu:ngpus=1
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR


# Pressure equilibration
pmemd.cuda \
    -O \
    -i equil.in \
    -o equil.out \
    -p initial_wat_ion.parm7 \
    -c heat.ncrst \
    -r equil.ncrst \
    -x equil.nc

# Production run
pmemd.cuda \
    -O \
    -i production.in \
    -o production.out \
    -p initial_wat_ion.parm7 \
    -c equil.ncrst \
    -r production.ncrst \
    -x production.nc

exit 0
