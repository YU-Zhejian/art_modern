#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=build
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

# pixi run -e prev art_profile_builder \
#     --i-file raw_data/ERR2888331.fastq \
#     --i-format FASTQ \
#     --read_len 400 \
#     --o-file1 DNBSeqG400L400R1.txt \
#     --parallel 4
# art-profile-fastqc --input DNBSeqG400L400R1.txt --output DNBSeqG400L400R1.png

pixi run -e prev art_profile_builder \
    --i-file raw_data/MGISEQ2000_PCR-free_NA12878_1_V100003043_L01_1.fq.gz \
    --i-format FASTQ \
    --read_len 150 \
    --o-file1 DNBSeqG400L150R1.txt \
    --parallel 4
art-profile-fastqc --input DNBSeqG400L150R1.txt --output DNBSeqG400L150R1.png
pixi run -e prev art_profile_builder \
    --i-file raw_data/MGISEQ2000_PCR-free_NA12878_1_V100003043_L01_2.fq.gz \
    --i-format FASTQ \
    --read_len 150 \
    --o-file1 DNBSeqG400L150R2.txt \
    --parallel 4
art-profile-fastqc --input DNBSeqG400L150R2.txt --output DNBSeqG400L150R2.png
echo "DONE"
