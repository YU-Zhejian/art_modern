#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem=2G
#SBATCH --time=4:00:00
#SBATCH --job-name=generate
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

mkdir -p generated

echo "Generating large contigs..."

pixi run -e testextreme \
    ./c/bin/generate_rand_large_contigs |
    pv >generated/large_contigs.fa

samtools faidx generated/large_contigs.fa
