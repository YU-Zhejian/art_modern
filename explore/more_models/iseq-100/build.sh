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

pixi run -e prev art_profile_builder \
    --i-file raw_data/SRR8791276_1.fastq \
    --i-format FASTQ \
    --read_len 150 \
    --o-file1 ISeq100L100R1.txt \
    --parallel 4

pixi run -e prev art_profile_builder \
    --i-file raw_data/SRR8791276_2.fastq \
    --i-format FASTQ \
    --read_len 150 \
    --o-file1 ISeq100L100R2.txt \
    --parallel 4
art-profile-fastqc --input ISeq100L100R1.txt --output ISeq100L100R1.png
art-profile-fastqc --input ISeq100L100R2.txt --output ISeq100L100R2.png
echo "DONE"
# | `ISeq100_150bp` | `ISeq100L100R1` | Y | 150 | 150 | iSeq100 | v1.4.0 | [SRR8791276](https://trace.ncbi.nlm.nih.gov/Traces/sra?run=SRR8791276), [DOI](https://doi.org/10.1371/journal.pmed.1002823) |
