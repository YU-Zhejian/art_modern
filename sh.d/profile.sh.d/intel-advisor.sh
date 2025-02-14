#!/usr/bin/env bash
# shellcheck disable=SC2317
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
PROFILE_DIR=opt/build_profile_intel_advisor
mkdir -p "${PROFILE_DIR}"

env -C "${PROFILE_DIR}" cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=ONEMKL \
    -G Ninja "$(pwd)"

env -C "${PROFILE_DIR}" ninja -j120

"${PROFILE_DIR}"/art_modern --version
ldd "${PROFILE_DIR}"/art_modern

advisor \
    --collect=roofline \
    --project-dir="${PROFILE_DIR}"/advi_results -- \
    "${PROFILE_DIR}"/art_modern \
    --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 data/Illumina_profiles/HiSeq2500L125R2.txt \
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
    --pe_frag_dist_mean 500 \
    --o-sam "${PROFILE_DIR}"/test_small_se_wgs_memory.bam \
    --o-sam-write_bam \
    --o-hl_sam "${PROFILE_DIR}"/test_small_se_wgs_memory.hl.sam \
    --o-fastq "${PROFILE_DIR}"/test_small_se_wgs_memory.fastq \
    --o-pwa "${PROFILE_DIR}"/test_small_se_wgs_memory.pwa

advisor-gui "${PROFILE_DIR}"/advi_results
