#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 1
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
# env -C raw_data pixi run -e moremodels axel \
#     https://ftp.cngb.org/pub/CNSA/data7/CNP0008318/CNS1487701/CNX1327094/CNR1513487/DL100002836_L02_UDB_1.fq.gz
env -C raw_data pixi run -e moremodels axel \
    https://ftp.cngb.org/pub/CNSA/data7/CNP0008318/CNS1487701/CNX1327094/CNR1513487/DL100002836_L02_UDB_2.fq.gz

echo "DONE"
