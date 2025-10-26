#!/usr/bin/env bash
CFLAGS=(-O3 -mtune=native -march=native -Wall -Wextra -Wpedantic)
CC="${CC:-cc}"
set -ue
SHDIR="$(readlink -f "$(dirname "${0}")")"
cd "${SHDIR}/../"
"${CC}" "${CFLAGS[@]}" -std=c11 \
    -o "${SHDIR}"/test-am-many-contigs.d/generate_many_contigs \
    "${SHDIR}"/test-am-many-contigs.d/generate_many_contigs.c

"${SHDIR}"/test-am-many-contigs.d/generate_many_contigs |
    opt/build_release_install/bin/art_modern \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file /dev/stdin \
        --i-type fasta \
        --read_len 125 \
        --mode template \
        --lc se \
        --i-parser stream \
        --i-fcov 2 \
        --parallel 0
# First, we will not attach any writer
exit

--o-hl_sam /dev/null \
    --o-hl_sam-num_threads 4 \
    --o-hl_sam-compress_level u \
    --o-hl_sam-write_bam \
    --o-fastq /dev/null
