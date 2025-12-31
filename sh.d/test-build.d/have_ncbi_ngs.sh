#!/usr/bin/env bash
set -uex
cd "$(readlink -f "$(dirname "$0")")"

latf-load --no-readnames -p ILLUMINA 1_1.fq 1_2.fq -o out.sra.d -L info -q PHRED_33
kar -f -c out.sra -d out.sra.d --md5
kar --long-list --test out.sra
fastq-dump ./out.sra -O out.sra.fastq.d --split-3

g++ \
    test_ncbi_ngs.cc \
    -lncbi-ngs-c++ \
    -lncbi-ngs \
    -lngs-c++ \
    -lncbi-vdb \
    -o test_ncbi_ngs
./test_ncbi_ngs ./out.sra
rm -fr test_ncbi_ngs out.sra out.sra.md5 out.sra.d out.sra.fastq.d
