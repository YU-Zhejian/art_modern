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
env -C raw_data axel \
    https://downloads.pacbcloud.com/public/onso/2023Q3/WGS/hg002_30x_WGS/Onso_hg002_PCR_free_WGS_OSQ_R1.fastq.gz
env -C raw_data axel \
    https://downloads.pacbcloud.com/public/onso/2023Q3/WGS/hg002_30x_WGS/Onso_hg002_PCR_free_WGS_OSQ_R2.fastq.gz

echo "DONE"
