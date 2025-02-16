#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set -ue
PROFILE_DIR=opt/build_release

# make release

"${PROFILE_DIR}"/art_modern --version
ldd "${PROFILE_DIR}"/art_modern

perf record -F 250 -a -g -- \
    "${PROFILE_DIR}"/art_modern \
    --builtin_qual_file HiSeq2500_125bp \
    --i-file data/raw_data/ce11_chr1.fa \
    --read_len 125 \
    --mode wgs \
    --lc pe \
    --i-parser memory \
    --i-fcov 10 \
    --parallel 0 \
    --ins_rate_1 0.1 \
    --del_rate_1 0.1 \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500
