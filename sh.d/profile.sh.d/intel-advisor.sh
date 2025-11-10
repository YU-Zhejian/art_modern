#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
PROFILE_DIR=opt/build_profile_intel_advisor
mkdir -p "${PROFILE_DIR}"

env -C "${PROFILE_DIR}" cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=ON \
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
    --i-file data/raw_data/ce11_chr1.fa \
    --mode wgs \
    --lc pe \
    --i-parser memory \
    --i-fcov 200 \
    --parallel 20 \
    --ins_rate_1 0.1 \
    --del_rate_1 0.1 \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500

advisor-gui "${PROFILE_DIR}"/advi_results
