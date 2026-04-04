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
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/MGISEQ/NA12878_1/MGISEQ2000_PCR-free_NA12878_1_V100003043_L01_1.fq.gz
env -C raw_data pixi run -e moremodels axel \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/MGISEQ/NA12878_1/MGISEQ2000_PCR-free_NA12878_1_V100003043_L01_2.fq.gz

pixi run -e moremodels fasterq-dump --progress -O raw_data/ ERR2888331

echo "DONE"
