#!/usr/bin/env bash
# This checks the bug that is found in 1.3.0 and 1.3.1.
set -ue
SHDIR="$(readlink -f "$(dirname "${0}")")"
cd "${SHDIR}/../"

export OPT_DIR="${PWD}/opt/test_until_failure"
export CC=icx
export CXX=icpx
export WORKDIR="${OPT_DIR}/workdir"

make release \
    CMAKE_FLAGS='-DUSE_RANDOM_GENERATOR=PCG -DUSE_MALLOC=JEMALLOC -DUSE_THREAD_PARALLEL=NOP'

NUM_TESTS=0
while true; do
    NUM_TESTS=$((NUM_TESTS + 1))
    echo "Starting test #${NUM_TESTS}..."
    rm -fr "${WORKDIR}"
    mkdir -p "${WORKDIR}"
   
        "${OPT_DIR}/build_release_install/bin/art_modern" \
        --i-file "${PWD}"/data/raw_data/ce11.mRNA_head.pbsim3.transcript \
        --i-type pbsim3_transcripts \
        --i-batch_size 100 \
        --mode template \
        --lc pe \
        --i-parser stream \
        --parallel 0 \
        --read_len_1 10 \
        --read_len_2 150 \
        --ins_rate_1 0.1 \
        --del_rate_1 0.1 \
        --o-hl_sam "${WORKDIR}"/test_small_pe_template_stream_pbsim3.hl.sam \
        --o-fastq "${WORKDIR}"/test_small_pe_template_stream_pbsim3.fq
    # Stop at 200th test to avoid infinite loop during demonstration
    if [ "${NUM_TESTS}" -ge 200 ]; then
        echo "Reached ${NUM_TESTS} tests without failure. Stopping."
        break
    fi
done
