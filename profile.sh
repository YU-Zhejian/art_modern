#!/usr/bin/env bash
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
# rm -fr build_profile
mkdir -p build_profile
cd build_profile
cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -G Ninja ..
ninja -j120
./art_modern --version
ldd ./art_modern

advisor \
    --collect=roofline \
    --project-dir=./advi_results -- \
    ./art_modern \
    --qual_file_1 ../art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 ../art/Illumina_profiles/HiSeq2500L125R2.txt \
    --i-file ../raw_data/ce11_chr1.fa \
    --read_len 125 \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 20 \
    --parallel 0 \
    --ins_rate_1 0.1 \
    --del_rate_1 0.1 \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500
advisor-gui ./advi_results
cd ..
#
#    --o-sam test_small_se_wgs_memory.sam \
#    --o-sam-write_bam \
