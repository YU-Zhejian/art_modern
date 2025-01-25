#!/usr/bin/env bash
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

# Build HTSLib
env -C src/htslib-1.21 \
    ./configure --prefix="$(pwd)"/opt \
    CC=icx \
    CFLAGS='-Ofast -mtune=native -march=native'
env -C src/htslib-1.21 make -j20
env -C src/htslib-1.21 make -j20 install

env -C src/gsl-2.8 \
    ./configure --prefix="$(pwd)"/opt \
    --enable-shared=yes \
    --enable-static=yes \
    CC=icx \
    CFLAGS='-Ofast -mtune=native -march=native'
env -C src/gsl-2.8 make -j20
env -C src/gsl-2.8 make -j20 install

# Build Original ART
icpx -Ofast -w -mtune=native -march=native \
    -lgsl -lgslcblas \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -o bin/art_original src/art_original/*.cpp

# Build wgsim
icpx -Ofast -w -mtune=native -march=native \
    -lhts -lz -lpthread -lm -lc \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -o bin/wgsim src/wgsim.c

# Build DWGSIM
icx -Ofast -w -mtune=native -march=native \
    -lhts -lm -lc \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -D_FILE_OFFSET_BITS=64 \
    -D_LARGEFILE64_SOURCE \
    -D_USE_KNETFILE \
    -DPACKAGE_VERSION='"0.1.15"' \
    -o bin/dwgsim src/dwgsim/*.c

icpx -Ofast -w -fopenmp -std=c++17 -march=native -mtune=native \
    -lz -lpthread \
    -DSFMT_MEXP=19937 -DHAVE_CONFIG_H -DPKGDATADIR='"/usr/local/share/pirs"' \
    -Isrc/pirs/SFMT-src-1.4 \
    -o bin/pirs \
    src/pirs/*.cpp src/pirs/SFMT-src-1.4/SFMT.c

# Build previous version of ART_MODERN
rm -fr src/am_prev_ver
git clone "$(git remote get-url origin)" src/am_prev_ver -b master
mkdir -p opt/art_modern_prev_ver_build/
env -C opt/art_modern_prev_ver_build/ cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_THREAD_PARALLEL=ASIO \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=ONEMKL \
    -DUSE_HTSLIB=hts \
    -DCMAKE_PREFIX_PATH="$(pwd)"/opt \
    -DCMAKE_INCLUDE_PATH="$(pwd)"/opt/include \
    -G Ninja "$(pwd)"/src/am_prev_ver
env -C opt/art_modern_prev_ver_build/ ninja
