#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --job-name=buildoncluster
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

# SUBMIT THIS JOB ON ROOT DIRECTORY OF THE REPO.

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

pwd
pixi run --frozen -e buildoncluster make release-mpi \
    CMAKE_FLAGS='-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx' \
    JOBS=16
