#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -N 1
# Default mem 4GiB
#SBATCH --time=3-00:00:00
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
    https://ftp.cngb.org/pub/CNSA/data2/CNP0000675/CNS0111621/CNX0095122/CNR0117350/T10_NA24695-2_1.fq.gz
env -C raw_data pixi run -e moremodels axel \
    https://ftp.cngb.org/pub/CNSA/data2/CNP0000675/CNS0111621/CNX0095122/CNR0117350/T10_NA24695-2_2.fq.gz

echo "DONE"
