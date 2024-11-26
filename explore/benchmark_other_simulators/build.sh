#!/usr/bin/env bash
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

# Build HTSLib
env -C src/htslib-1.21 \
    ./configure --prefix="$(pwd)"/opt \
    CC=icx \
    CFLAGS='-O3 -mtune=native'
env -C src/htslib-1.21 make -j20
env -C src/htslib-1.21 make -j20 install

# Build art_modern
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

# Build Original ART
icpx -O3 -w -mtune=native \
    -lgsl \
    -o bin/art_original src/art_original/*.cpp

# Build wgsim
icpx -O3 -w -mtune=native \
    -lhts -lz -lpthread -lm -lc \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -o bin/wgsim src/wgsim.c
