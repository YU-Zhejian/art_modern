#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
AMD_COMPILER_PATH="${AMD_COMPILER_PATH:-/opt/AMD/aocc-compiler-5.0.0/}"
AMD_UPROF_PATH="${AMD_COMPILER_PATH:-/opt/AMDuProf_5.0-1479/}"
AMD_AOCL_PATH="${AMD_COMPILER_PATH:-/opt/AMD/aocl/aocl-linux-aocc-5.0.0/aocc/}"

if [ ! -d "${AMD_COMPILER_PATH}" ]; then
    echo "AMD compiler not found at ${AMD_COMPILER_PATH}"
    exit 1
else
    echo "AMD compiler found at ${AMD_COMPILER_PATH}"
fi

if [ ! -d "${AMD_UPROF_PATH}" ]; then
    echo "AMD uprof not found at ${AMD_UPROF_PATH}"
    exit 1
else
    echo "AMD uprof found at ${AMD_UPROF_PATH}"
fi

if [ ! -d "${AMD_AOCL_PATH}" ]; then
    echo "AMD AOCL not found at ${AMD_AOCL_PATH}"
    exit 1
else
    echo "AMD AOCL found at ${AMD_AOCL_PATH}"
fi

set +ue
. "${AMD_COMPILER_PATH}"/setenv_AOCC.sh
. "${AMD_AOCL_PATH}"/amd-libs.cfg
set -ue

PROFILE_DIR=opt/build_profile_amd_uprof
mkdir -p "${PROFILE_DIR}"

env -C "${PROFILE_DIR}" cmake \
    -DCMAKE_C_COMPILER="${AMD_COMPILER_PATH}"/bin/clang \
    -DCMAKE_CXX_COMPILER="${AMD_COMPILER_PATH}"/bin/clang++ \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=STL \
    -G Ninja "$(pwd)"

env -C "${PROFILE_DIR}" ninja -j120

"${PROFILE_DIR}"/art_modern --version
ldd "${PROFILE_DIR}"/art_modern
rm -fr "${PROFILE_DIR}"/uprof_results
"${AMD_UPROF_PATH}"/bin/AMDuProfCLI \
    collect \
    --config hotspots \
    --output-dir "${PROFILE_DIR}"/uprof_results \
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
"${AMD_UPROF_PATH}"/bin/AMDuProf \
    --session "${PROFILE_DIR}"/uprof_results/* \
    --bin-path "${PROFILE_DIR}" --src-path "$(pwd)"
