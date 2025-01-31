#!/usr/bin/env bash
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
rm -fr opt/art_modern_build/
mkdir -p opt/art_modern_build/
env -C opt/art_modern_build/ cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_THREAD_PARALLEL=BS \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=ONEMKL \
    -DUSE_HTSLIB=hts \
    -DCMAKE_PREFIX_PATH="$(pwd)"/opt \
    -DCMAKE_INCLUDE_PATH="$(pwd)"/opt/include \
    -G Ninja "$(pwd)"/../../
env -C opt/art_modern_build/ ninja

rm -fr opt/art_modern_gcc_build/
mkdir -p opt/art_modern_gcc_build/
env -C opt/art_modern_gcc_build/ cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_THREAD_PARALLEL=BS \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=STL \
    -G Ninja "$(pwd)"/../../
env -C opt/art_modern_gcc_build/ ninja

rm -fr opt/art_modern_pcg_build/
mkdir -p opt/art_modern_pcg_build/
env -C opt/art_modern_pcg_build/ cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_THREAD_PARALLEL=BS \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=PCG \
    -DUSE_HTSLIB=hts \
    -DCMAKE_PREFIX_PATH="$(pwd)"/opt \
    -DCMAKE_INCLUDE_PATH="$(pwd)"/opt/include \
    -G Ninja "$(pwd)"/../../
env -C opt/art_modern_pcg_build/ ninja
