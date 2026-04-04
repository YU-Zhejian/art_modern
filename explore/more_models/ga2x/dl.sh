#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem=4G
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
pixi run -e moremodels fasterq-dump --split-files --progress --outdir raw_data/ SRR933371
pixi run -e moremodels fasterq-dump --split-files --progress --outdir raw_data/ ERR038223

echo "DONE"
