#!/usr/bin/env bash
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
mkdir -p opt/art_modern_build/
env -C opt/art_modern_build/ cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=ONEMKL \
    -DUSE_HTSLIB=hts \
    -DCMAKE_PREFIX_PATH="$(pwd)"/opt \
    -DCMAKE_INCLUDE_PATH="$(pwd)"/opt/include \
    -G Ninja "$(pwd)"/../../
env -C opt/art_modern_build/ ninja
