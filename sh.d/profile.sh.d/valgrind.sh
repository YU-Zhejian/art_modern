#!/usr/bin/env bash
set -ue
PROFILE_DIR=opt/build_profile_valgrind
mkdir -p "${PROFILE_DIR}"

env -C "${PROFILE_DIR}" cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -Wdev -Wdeprecated --warn-uninitialized \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_RANDOM_GENERATOR=PCG \
    -DUSE_MALLOC=JEMALLOC \
    -G Ninja "$(pwd)"

env -C "${PROFILE_DIR}" ninja -j120

VALGRIND_LOG="${PROFILE_DIR}/valgrind.log"

valgrind \
    --tool=memcheck \
    --leak-check=full \
    --show-leak-kinds=all \
    --track-origins=yes \
    --errors-for-leak-kinds=definite,possible \
    --num-callers=40 \
    --malloc-fill=0xAA \
    --free-fill=0xDA \
    --gen-suppressions=all \
    --log-file="${VALGRIND_LOG}" \
    "${PROFILE_DIR}"/art_modern \
    --i-file data/raw_data/lambda_phage.fa \
    --mode wgs \
    --lc pe \
    --i-parser memory \
    --i-fcov 4 \
    --parallel 2 \
    --ins_rate_1 0.1 \
    --del_rate_1 0.1 \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500 \
    --o-sam "${PROFILE_DIR}"/test_small_se_wgs_memory.bam \
    --o-sam-write_bam \
    --o-hl_sam "${PROFILE_DIR}"/test_small_se_wgs_memory.hl.sam \
    --o-fastq "${PROFILE_DIR}"/test_small_se_wgs_memory.fastq \
    --o-pwa "${PROFILE_DIR}"/test_small_se_wgs_memory.pwa
