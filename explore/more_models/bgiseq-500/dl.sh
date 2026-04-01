#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=raw_data-dl
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

# prefetch --progress -O raw_data SRR36802437 --max-size 0
fasterq-dump --progress -O raw_data SRR36802437 --split-files
echo "DONE"
