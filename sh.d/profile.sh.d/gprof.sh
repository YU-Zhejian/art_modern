#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091

# See: <https://cygwin.com/cygwin-ug-net/gprof.html>
set -uex
PROFILE_DIR="$(pwd)/opt/build_profile_gprof"
mkdir -p "${PROFILE_DIR}"

. sh.d/libwindows.sh

env -C "${PROFILE_DIR}" \
    CFLAGS='-pg' \
    CXXFLAGS='-pg' \
    LDFLAGS='-pg -lgmon' \
    cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=PCG \
    ${CMAKE_FLAGS:-} \
    -G Ninja "$(pwd)"

env -C "${PROFILE_DIR}" ninja -j120 -v

if is_windows; then
    # Move DLLs
    cp "${PROFILE_DIR}"/deps/*/cyg*.dll "${PROFILE_DIR}"
fi

"${PROFILE_DIR}"/art_modern --version
ldd "${PROFILE_DIR}"/art_modern
nm "${PROFILE_DIR}"/art_modern | grep gmon # Check for gmon symbols

export GMON_OUT_PREFIX="${PROFILE_DIR}"/gmon
"${PROFILE_DIR}"/art_modern \
    --i-file data/raw_data/ce11_chr1.fa \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 2 \
    --parallel 20
gprof "${PROFILE_DIR}"/art_modern "${PROFILE_DIR}"/gmon.out* > "${PROFILE_DIR}"/gprof_report.txt
