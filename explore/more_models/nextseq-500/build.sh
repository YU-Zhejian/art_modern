#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --job-name=build
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

pixi run -e prev art_profile_builder \
    --i-file raw_data/SRR9733244_1.fastq \
    --i-format FASTQ \
    --read_len 150 \
    --o-file1 NextSeq500L150R1.txt \
    --parallel 4

pixi run -e prev art_profile_builder \
    --i-file raw_data/SRR9733244_2.fastq \
    --i-format FASTQ \
    --read_len 150 \
    --o-file1 NextSeq500L150R2.txt \
    --parallel 4
art-profile-fastqc --input NextSeq500L150R1.txt --output NextSeq500L150R1.png
art-profile-fastqc --input NextSeq500L150R2.txt --output NextSeq500L150R2.png

echo "DONE"
