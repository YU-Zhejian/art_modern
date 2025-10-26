#!/usr/bin/env bash
CFLAGS=(-O3 -mtune=native -march=native -Wall -Wextra -Wpedantic)
CC="${CC:-cc}"
set -ue
SHDIR="$(readlink -f "$(dirname "${0}")")"
cd "${SHDIR}/../"
make -C "${SHDIR}"/test-am-large-contigs.d/


    opt/build_release_install/bin/art_modern \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file "${SHDIR}"/test-am-large-contigs.d/large_contigs.fa.gz \
        --i-type fasta \
        --read_len 125 \
        --mode wgs \
        --lc se \
        --i-parser htslib \
        --i-fcov 20 \
        --parallel 0 \
        --o-sam /dev/null \
        --o-sam-num_threads 4 \
        --o-sam-compress_level u \
        --o-sam-write_bam \
        --o-fastq /dev/null
