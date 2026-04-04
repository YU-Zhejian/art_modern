#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=generate
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

NJOBS=32
mkdir -p generated

echo "Generating large contigs..."

pixi run -e testextreme \
    ./c/bin/generate_rand_large_contigs > generated/large_contigs.fa

samtools faidx generated/large_contigs.fa

bgzip --index --compress-level 9 --keep -@"${NJOBS}" generated/large_contigs.fa
samtools faidx generated/large_contigs.fa.gz
