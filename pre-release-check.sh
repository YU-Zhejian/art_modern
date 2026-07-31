#!/usr/bin/env bash
# shellcheck disable=SC1091
# SC2317: setvars.sh ends in `return`, which `shellcheck -x` models as terminating this script, so
# every command below the `.` reads as dead code. A `return` in a sourced file only ends the `.`.
# shellcheck disable=SC2317
set -ue

# WARNING: Have hard-coded paths.
# To tell whether the build is successful, see whether there are leftover test-build-* directories in the current working directory. If there are, the build failed.

# Running overnight is highly recommended, as this script will take a long time to run.

# This script sits in the project root. cd there so that `make`, the log-dir cleanup and
# CMAKE_TOOLCHAIN_FILE below do not depend on where the script was invoked from.
PROJDIR="$(readlink -f "$(dirname "${0}")")"
cd "${PROJDIR}"

# Everything printed from here on is also written to this file. A test-build.py that aborts before
# creating its log dir leaves no test-build-*.log.d behind, so without this its error is lost and the
# run looks clean.
PRC_LOG="${PROJDIR}/pre-release-check-$(date +%Y%m%d-%H%M%S).log"
exec > >(tee -a "${PRC_LOG}") 2>&1
echo "Logging this run to ${PRC_LOG}"

# Run GCC: Run GCC-compatible links. ICX/ICPX belong here too; only Clang links LLVM-style, which is
# why run_llvm needs its own Boost and {fmt} built against libc++.
function run_gcc() {
    LD_LIBRARY_PATH="${HOME}/opt/sra-tools-3.3.0/lib64/:${HOME}/opt/ncbi-vdb-3.3.0/lib64/:${LD_LIBRARY_PATH:-}" \
        LD_RUN_PATH="${HOME}/opt/sra-tools-3.3.0/lib64/:${HOME}/opt/ncbi-vdb-3.3.0/lib64/:${LD_RUN_PATH:-}" \
        LIBRARY_PATH="${HOME}/opt/sra-tools-3.3.0/lib64/:${HOME}/opt/ncbi-vdb-3.3.0/lib64/:${LIBRARY_PATH:-}" \
        CMAKE_LIBRARY_PATH="${HOME}/opt/sra-tools-3.3.0/lib64/:${HOME}/opt/ncbi-vdb-3.3.0/lib64/:${CMAKE_LIBRARY_PATH:-}" \
        PATH="${HOME}/opt/sra-tools-3.3.0/bin/:${HOME}/opt/ncbi-vdb-3.3.0/bin/:${PATH:-}" \
        C_INCLUDE_PATH="${HOME}/opt/sra-tools-3.3.0/include/:${HOME}/opt/ncbi-vdb-3.3.0/include/:${C_INCLUDE_PATH:-}" \
        CPLUS_INCLUDE_PATH="${HOME}/opt/sra-tools-3.3.0/include/:${HOME}/opt/ncbi-vdb-3.3.0/include/:${CPLUS_INCLUDE_PATH:-}" \
        "${@}"
}

function run_llvm() {
    # SRA-Tools and NCBI-VDB are not yet built with LLVM, so we don't set those paths here.
    LD_LIBRARY_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${HOME}/opt/fmt-12.0.0-clang/lib/:${LD_LIBRARY_PATH:-}" \
        LD_RUN_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${HOME}/opt/fmt-12.0.0-clang/lib/:${LD_RUN_PATH:-}" \
        CMAKE_PREFIX_PATH="${HOME}/opt/fmt-12.0.0-clang/lib/cmake/:${CMAKE_PREFIX_PATH:-}" \
        PKG_CONFIG_PATH="${HOME}/opt/fmt-12.0.0-clang/lib/pkgconfig/:${PKG_CONFIG_PATH:-}" \
        CMAKE_TOOLCHAIN_FILE="${PROJDIR}/sh.d/toolchain/host-llvm/llvm-toolchain.cmake" \
        ASSERT_USING_LLVM_CXXSTDLIB=1 \
        "${@}"
}

# Checked before anything else: the Intel stages come last, and without this the whole overnight run
# would be spent before finding out.
if [ ! -f "/opt/intel/oneapi/setvars.sh" ]; then
    echo "Intel OneAPI not detected."
    exit 1
fi

make clean
# `make clean` does not cover these, and a directory left by an earlier run would be misread as a
# failure of this one.
rm -fr test-build-*.log.d

# Modify this line for whether to testsmall.
TB_TARGET="testbuild-small" #

run_gcc make "${TB_TARGET}" || true
run_gcc make "${TB_TARGET}"-mpi || true

run_llvm make "${TB_TARGET}" \
    CMAKE_FLAGS="-DBoost_DIR=${HOME}/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/" || true

run_llvm make "${TB_TARGET}"-mpi \
    CMAKE_FLAGS="-DBoost_DIR=${HOME}/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/" || true

set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
run_gcc make "${TB_TARGET}" \
    CMAKE_FLAGS='-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx' || true
run_gcc make "${TB_TARGET}"-mpi \
    CMAKE_FLAGS='-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx' || true
