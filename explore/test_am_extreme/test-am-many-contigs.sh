#!/usr/bin/env bash
CFLAGS=(-O3 -mtune=native -march=native -Wall -Wextra -Wpedantic)
CC="${CC:-cc}"
set -ue
SHDIR="$(readlink -f "$(dirname "${0}")")"
cd "${SHDIR}/../"
make -C "${SHDIR}"/test-am-many-contigs.d/

"${SHDIR}"/test-am-many-contigs.d/generate_many_contigs |
    opt/build_release_install/bin/art_modern \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file /dev/stdin \
        --i-type fasta \
        --mode template \
        --lc se \
        --i-parser stream \
        --i-fcov 2 \
        --i-batch_size 1048576 \
        --parallel 0 \
        --o-hl_sam /dev/null \
        --o-hl_sam-num_threads 4 \
        --o-hl_sam-compress_level u \
        --o-hl_sam-write_bam \
        --o-fastq /dev/null \
        --reporting_interval-job_executor 10 \
        --reporting_interval-job_pool 50
# This generates 10G reads.
