#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091

# Running overnight is highly recommended, as this script will take a long time to run.

function run_gcc() {
    LD_LIBRARY_PATH="${HOME}/opt/sra-tools-3.3.0/lib64/:${HOME}/opt/ncbi-vdb-3.3.0/lib64/:${LD_LIBRARY_PATH:-}" \
        LD_RUN_PATH="${HOME}/opt/sra-tools-3.3.0/lib64/:${HOME}/opt/ncbi-vdb-3.3.0/lib64/:${LD_RUN_PATH:-}" \
        LIBRARY_PATH="${HOME}/opt/sra-tools-3.3.0/lib64/:${HOME}/opt/ncbi-vdb-3.3.0/lib64/:${LIBRARY_PATH:-}" \
        PATH="${HOME}/opt/sra-tools-3.3.0/bin/:${HOME}/opt/ncbi-vdb-3.3.0/bin/:${PATH:-}" \
        C_INCLUDE_PATH="${HOME}/opt/sra-tools-3.3.0/include/:${HOME}/opt/ncbi-vdb-3.3.0/include/:${C_INCLUDE_PATH:-}" \
        CPLUS_INCLUDE_PATH="${HOME}/opt/sra-tools-3.3.0/include/:${HOME}/opt/ncbi-vdb-3.3.0/include/:${CPLUS_INCLUDE_PATH:-}" \
        "${@}"
}

function run_llvm() {
    # SRA-Tools and NCBI-VDB are not yet built with LLVM, so we don't set those paths here.
    LD_LIBRARY_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${HOME}/opt/fmt-12.0.0-clang/lib/:${LD_LIBRARY_PATH:-}" \
        LD_RUN_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${HOME}/opt/fmt-12.0.0-clang/lib/:${LD_RUN_PATH:-}" \
        PKG_CONFIG_PATH="${HOME}/opt/fmt-12.0.0-clang/lib/pkgconfig/:${PKG_CONFIG_PATH:-}" \
        CMAKE_TOOLCHAIN_FILE="$(pwd)/sh.d/toolchain/host-llvm/llvm-toolchain.cmake" \
        ASSERT_USING_LLVM_CXXSTDLIB=1 \
        "${@}"
}

make clean

# run_gcc make testbuild-small || true
# run_gcc make testbuild-small-mpi || true

run_llvm make testbuild-small \
    CMAKE_FLAGS="-DBoost_DIR=${HOME}/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/" || true

run_llvm make testbuild-small-mpi \
    CMAKE_FLAGS="-DBoost_DIR=${HOME}/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/" || true

. /opt/intel/oneapi/setvars.sh
run_gcc make testbuild-small \
    CMAKE_FLAGS='-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx' || true
run_gcc make testbuild-small-mpi \
    CMAKE_FLAGS='-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx' || true
