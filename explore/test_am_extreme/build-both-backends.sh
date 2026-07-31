#!/usr/bin/bash
set -ueo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(cd "${HERE}/../.." && pwd)"
OPT_DIR="${HERE}/opt"
BUILD_TYPE="${BUILD_TYPE:-Debug}"
JOBS="${JOBS:-$(nproc)}"

build_backend() {
    local name="$1" backend="$2"
    local build_dir="${OPT_DIR}/build_am_${name}"
    local install_dir="${OPT_DIR}/am_${name}"

    rm -rf "${build_dir}" "${install_dir}"
    mkdir -p "${build_dir}"
    env -C "${build_dir}" cmake \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DUSE_THREAD_PARALLEL="${backend}" \
        -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
        -DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
        -DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
        -DCMAKE_INSTALL_PREFIX="${install_dir}" \
        "${SRC_DIR}"
    cmake --build "${build_dir}" -j"${JOBS}"
    cmake --install "${build_dir}"
}

assert_backend() {
    local name="$1" expected="$2"
    local exe="${OPT_DIR}/am_${name}/bin/art_modern"
    local version
    version="$("${exe}" --version 2>&1)"
    grep -qi "BS::thread_pool: ${expected}" <<<"${version}" || {
        echo "am_${name}: wrong backend" >&2
        exit 1
    }
    echo "am_${name}: OK ${exe}"
}

build_backend asio ASIO
build_backend bs BS
assert_backend asio "not used"
assert_backend bs "[0-9]"
