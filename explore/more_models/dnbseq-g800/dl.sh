#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -N 1
# Default mem 4GiB
#SBATCH --time=24:00:00
#SBATCH --job-name=raw_data-dl
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

mkdir -p raw_data
env -C raw_data pixi run -e moremodels axel \
    https://demodata.completegenomics.mgiamericas.com/Demo_Data/G800/G800_WGS_PE150_Ecoli_Read1.fq.gz
env -C raw_data pixi run -e moremodels axel \
    https://demodata.completegenomics.mgiamericas.com/Demo_Data/G800/G800_WGS_PE150_Ecoli_Read2.fq.gz

env -C raw_data pixi run -e moremodels axel \
    https://demodata.completegenomics.mgiamericas.com/Demo_Data/G800/G800_WGS_SE600_HG001_Read1.fq.gz
echo "DONE"
