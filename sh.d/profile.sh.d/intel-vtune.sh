#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
PROFILE_DIR=opt/build_profile_intel_vtune
mkdir -p "${PROFILE_DIR}"
env -C "${PROFILE_DIR}" cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_THREAD_PARALLEL=BS \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=ONEMKL \
    -G Ninja "$(pwd)"

env -C "${PROFILE_DIR}" ninja -j120

"${PROFILE_DIR}"/art_modern --version
ldd "${PROFILE_DIR}"/art_modern

for collect in hotspots threading memory-consumption; do # hpc-performance memory-access io
    rm -fr "${PROFILE_DIR}"/vtune-"${collect}"
    vtune \
        -collect="${collect}" \
        -source-search-dir=".." \
        -result-dir="${PROFILE_DIR}"/vtune-"${collect}" -- \
        "${PROFILE_DIR}"/art_modern \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file data/raw_data/ce11_chr1.fa \
        --read_len 125 \
        --mode wgs \
        --lc pe \
        --i-parser memory \
        --i-fcov 100 \
        --parallel 20 \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 500 \
        --o-sam /dev/null \
        --o-sam-write_bam \
        --o-sam-num_threads 20 \
        --o-hl_sam /dev/null \
        --o-fastq /dev/null \
        --o-pwa /dev/null
    vtune-gui "${PROFILE_DIR}"/vtune-"${collect}"
done
