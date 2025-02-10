#!/usr/bin/env bash
# shellcheck disable=SC1091
# shellcheck disable=SC2317

set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

# Build previous version of ART_MODERN
rm -fr src/am_prev_ver
git clone "$(git remote get-url origin)" src/am_prev_ver
env -C src/am_prev_ver git checkout 1.1.1
env -C src/am_prev_ver git apply "$(pwd)/patches/1.1.1.patch"
mkdir -p opt/art_modern_prev_ver_build/
env -C opt/art_modern_prev_ver_build/ cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_THREAD_PARALLEL=BS \
    -DUSE_BTREE_MAP=OFF \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=ONEMKL \
    -DUSE_HTSLIB=hts \
    -DCMAKE_PREFIX_PATH="$(pwd)"/opt \
    -DCMAKE_INCLUDE_PATH="$(pwd)"/opt/include \
    -G Ninja "$(pwd)"/src/am_prev_ver
env -C opt/art_modern_prev_ver_build/ ninja
