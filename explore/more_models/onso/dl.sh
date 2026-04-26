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
    https://prod-dcd-datasets-public-files-eu-west-1.s3.eu-west-1.amazonaws.com/81c1b971-ec7e-4097-9ef7-dbc61d5c5b9f -o WS2442508A01_100_R1.fastq.gz
env -C raw_data pixi run -e moremodels axel \
    https://prod-dcd-datasets-public-files-eu-west-1.s3.eu-west-1.amazonaws.com/00079839-a531-4086-82b8-8e0ba3c6b9f0 -o WS2442508A01_100_R2.fastq.gz

echo "DONE"
