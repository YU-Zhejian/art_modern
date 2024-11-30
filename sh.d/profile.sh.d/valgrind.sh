#!/usr/bin/env bash
set -ue
PROFILE_DIR=opt/build_profile_valgrind
mkdir -p "${PROFILE_DIR}"

env -C "${PROFILE_DIR}" cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCEU_CM_SHOULD_USE_NATIVE=OFF \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_RANDOM_GENERATOR=STL \
    -G Ninja "$(pwd)"

env -C "${PROFILE_DIR}" ninja -j120

"${PROFILE_DIR}"/art_modern --version
ldd "${PROFILE_DIR}"/art_modern

valgrind \
    --tool=callgrind \
    --dump-instr=yes \
    --simulate-cache=yes \
    --collect-jumps=yes \
    --callgrind-out-file="${PROFILE_DIR}"/callgrind.out \
    "${PROFILE_DIR}"/art_modern \
    --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 data/Illumina_profiles/HiSeq2500L125R2.txt \
    --i-file data/raw_data/ce11_chr1.fa \
    --read_len 125 \
    --mode wgs \
    --lc pe \
    --i-parser memory \
    --i-fcov 4 \
    --parallel 0 \
    --ins_rate_1 0.1 \
    --del_rate_1 0.1 \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500 \
    --o-sam "${PROFILE_DIR}"/test_small_se_wgs_memory.bam \
    --o-sam-write_bam \
    --o-hl_sam "${PROFILE_DIR}"/test_small_se_wgs_memory.hl.sam \
    --o-fastq "${PROFILE_DIR}"/test_small_se_wgs_memory.fastq \
    --o-pwa "${PROFILE_DIR}"/test_small_se_wgs_memory.pwa
