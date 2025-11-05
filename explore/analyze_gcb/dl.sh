#!/usr/bin/env bash
set -ue

mkdir -p raw_data

# snRNA-seq of C. elegans: Adult Male Nuclei (SRR34148065)
#
for i in {1..2}; do
    curl \
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/065/SRR34148065/SRR34148065_"${i}".fastq.gz |
        head -n 1000000 >raw_data/SRR34148065_"${i}".fastq.gz
done
