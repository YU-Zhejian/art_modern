#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem=2GB
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
    https://ftp.cngb.org/pub/CNSA/data7/CNP0008217/CNS1414621/CNX1259793/CNR1433837/ML150002521_L01_UDB-386_1.fq.gz
env -C raw_data pixi run -e moremodels axel \
    https://ftp.cngb.org/pub/CNSA/data7/CNP0008217/CNS1414621/CNX1259793/CNR1433837/ML150002521_L01_UDB-386_2.fq.gz

echo "DONE"
