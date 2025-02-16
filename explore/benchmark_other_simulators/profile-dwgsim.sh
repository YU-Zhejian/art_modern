#!/usr/bin/env bash
# shellcheck disable=SC1091
# shellcheck disable=SC2317

set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
collect=hotspot
icx -Ofast -w -mtune=native -g \
    -lhts -lz -lm -lc \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -D_FILE_OFFSET_BITS=64 \
    -D_LARGEFILE64_SOURCE \
    -D_USE_KNETFILE \
    -DPACKAGE_VERSION='"0.1.15"' \
    -o bin/dwgsim-relwithdbginfo src/dwgsim/*.c
mkdir -p data/vtune-dwgsim-"${collect}"
vtune \
    -collect="${collect}" \
    -source-search-dir="." \
    -result-dir=data/vtune-dwgsim-"${collect}" -- \
    bin/dwgsim-relwithdbginfo \
    -1 150 -2 150 -C 10 -d 300 -s 20 -o 2 -r 0 -y 0 \
    data/ce11.fa \
    /tmp/ce11_dwgsim_profile
vtune-gui data/vtune-dwgsim-"${collect}"
rm -fr /tmp/ce11_dwgsim
