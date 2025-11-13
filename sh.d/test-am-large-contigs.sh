#!/usr/bin/env bash
set -ue
export CFLAGS=(-O3 -mtune=native -march=native -Wall -Wextra -Wpedantic)
export CC="${CC:-gcc}"
SHDIR="$(readlink -f "$(dirname "${0}")")"
cd "${SHDIR}/../"
make -C "${SHDIR}"/test-am-large-contigs.d/

# Only test whether the SAM output can be generated from such large contigs
/bin/time -a -o time.tsv -f '\t%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c' \
    opt/build_release_install/bin/art_modern \
    --builtin_qual_file HiSeq2500_125bp \
    --i-file "${SHDIR}"/test-am-large-contigs.d/large_contigs.fa.gz \
    --i-type fasta \
    --i-parser htslib \
    --i-fcov 20 \
    --parallel 0 \
    --o-sam /dev/null \
    --o-sam-num_threads 4 \
    --o-sam-write_bam \
    --o-fastq /dev/null
