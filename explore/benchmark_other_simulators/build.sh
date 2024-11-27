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
