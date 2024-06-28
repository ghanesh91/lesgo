#!/bin/bash
#PBS -A ULHI0003
#PBS -l walltime=12:00:00
#PBS -q main
#PBS -j oe
#PBS -m abe
#PBS -M naras062@umn.edu  
#PBS -l select=1:ncpus=12:mpiprocs=12

### Set TMPDIR as recommended
setenv TMPDIR /glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

module --force purge
module load ncarenv-basic/23.06
module load intel-oneapi/2023.0.0
module load intel-mpi/2021.8.0
module load hdf5-mpi/1.12.2
module load fftw-mpi/3.3.10
module load cmake

### Run the executable
mpiexec -np 12 ./lesgo-mpi >  lesgo_test.log
