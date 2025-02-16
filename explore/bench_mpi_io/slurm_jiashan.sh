#!/usr/bin/env bash
#SBATCH --job-name=test_mpi_io
#SBATCH --partition=cpu
#SBATCH --nodes=20
#SBATCH --mem=2GB
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=72:00:00

module () 
{ 
    eval `/usr/bin/modulecmd bash $*`
}


module load mpi/openmpi/4.1.4-gcc

set -ue
mpirun -np 20 opt/build/bench_one_io_per_job
