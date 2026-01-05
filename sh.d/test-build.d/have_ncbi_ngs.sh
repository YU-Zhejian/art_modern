#!/usr/bin/env bash
set -uex
cd "$(readlink -f "$(dirname "$0")")"
CXX="${CXX:-g++}"
DEFAULT_CXXFLAGS="-Wall -Wpedantic -std=c++11 -Og -g"
CXXFLAGS="${CXXFLAGS:-$DEFAULT_CXXFLAGS}"

if type latf-load >/dev/null 2>&1 && type kar >/dev/null 2>&1; then
    latf-load \
        --no-readnames \
        -p ILLUMINA \
        -o out.sra.d \
        -L info \
        -q PHRED_33 \
        1_1.fq 1_2.fq
    kar -f -c out.sra -d out.sra.d
else
    echo "latf-load or kar not found, using base64 version of out.sra"
    base64 -di out.sra.xz.base64 | xz -d >out.sra
fi

"${CXX}" \
    ${CXXFLAGS} \
    test_ncbi_ngs.cc \
    -lncbi-ngs \
    -lncbi-vdb \
    -o test_ncbi_ngs
./test_ncbi_ngs ./out.sra
rm -fr test_ncbi_ngs out.sra out.sra.d out.sra.fastq.d
