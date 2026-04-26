#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -N 1
#SBATCH --mem=32G
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
    --i-file raw_data/WS2442508A01_100_R1.fastq.gz \
    --i-format FASTQ \
    --read_len 150 \
    --o-file1 PacbOnsoL150R1.txt \
    --parallel 4
pixi run -e prev art_profile_builder \
    --i-file raw_data/WS2442508A01_100_R2.fastq.gz \
    --i-format FASTQ \
    --read_len 150 \
    --o-file1 PacbOnsoL150R2.txt \
    --parallel 4
art-profile-fastqc --input PacbOnsoL150R1.txt --output PacbOnsoL150R1.png
art-profile-fastqc --input PacbOnsoL150R2.txt --output PacbOnsoL150R2.png
echo "DONE"
