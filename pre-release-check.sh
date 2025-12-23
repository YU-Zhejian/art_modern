#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091

make clean
make testbuild-small || true
make testbuild-small-mpi || true

LD_LIBRARY_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${LD_LIBRARY_PATH:-}" \
    LD_RUN_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${LD_RUN_PATH:-}"\
    PKG_CONFIG_PATH="${HOME}/opt/fmt-12.0.0-clang/lib/pkgconfig/:${PKG_CONFIG_PATH:-}" \
    CMAKE_TOOLCHAIN_FILE="$(pwd)/sh.d/toolchain/host-llvm/llvm-toolchain.cmake" \
    make testbuild-small \
    CMAKE_FLAGS="-DBoost_DIR=${HOME}/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/" || true

LD_LIBRARY_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${LD_LIBRARY_PATH:-}" \
    LD_RUN_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${LD_RUN_PATH:-}"\
    PKG_CONFIG_PATH="${HOME}/opt/fmt-12.0.0-clang/lib/pkgconfig/:${PKG_CONFIG_PATH:-}" \
    CMAKE_TOOLCHAIN_FILE="$(pwd)/sh.d/toolchain/host-llvm/llvm-toolchain.cmake" \
    make testbuild-small-mpi \
    CMAKE_FLAGS="-DBoost_DIR=${HOME}/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/" || true

. /opt/intel/oneapi/setvars.sh
make testbuild-small \
    CMAKE_FLAGS='-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx' || true
make testbuild-small-mpi \
    CMAKE_FLAGS='-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx' || true


