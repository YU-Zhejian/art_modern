#!/usr/bin/env bash
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

# Build HTSLib
env -C src/htslib-1.21 \
    ./configure --prefix="$(pwd)"/opt \
    CC=icx \
    CFLAGS='-Ofast -mtune=native'
env -C src/htslib-1.21 make -j20
env -C src/htslib-1.21 make -j20 install

# Build Original ART
icpx -Ofast -w -mtune=native \
    -lgsl \
    -o bin/art_original src/art_original/*.cpp

# Build wgsim
icpx -Ofast -w -mtune=native \
    -lhts -lz -lpthread -lm -lc \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -o bin/wgsim src/wgsim.c

# Build DWGSIM
icx -Ofast -w -mtune=native \
    -lhts -lm -lc \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -D_FILE_OFFSET_BITS=64 \
    -D_LARGEFILE64_SOURCE \
    -D_USE_KNETFILE \
    -DPACKAGE_VERSION='"0.1.15"' \
    -o bin/dwgsim src/dwgsim/*.c
