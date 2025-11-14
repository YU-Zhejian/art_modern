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

item=roofline
advisor \
    --collect="${item}" \
    --enable-cache-simulation \
    --select 'has-source' \
    --project-dir="${PROFILE_DIR}"/advi_results-"${item}" -- \
    "${PROFILE_DIR}"/art_modern \
    --i-file data/raw_data/ce11_chr1.fa \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 2 \
    --parallel 20

advisor-gui "${PROFILE_DIR}"/advi_results-"${item}"
