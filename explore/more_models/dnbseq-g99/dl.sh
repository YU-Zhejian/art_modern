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
    https://ftp.cngb.org/pub/CNSA/data4/CNP0006051/CNS1176430/CNX1113182/CNR1261865/FT100044993_L01_101-101_1.fq.gz
env -C raw_data pixi run -e moremodels axel \
    https://ftp.cngb.org/pub/CNSA/data4/CNP0006051/CNS1176430/CNX1113182/CNR1261865/FT100044993_L01_101-101_2.fq.gz

echo "DONE"
