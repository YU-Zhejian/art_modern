#!/usr/bin/env bash
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
# rm -fr build_profile
mkdir -p build_profile
cd build_profile
cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
ninja -j120
./art_modern --version
ldd ./art_modern
advisor --report=roofline --project-dir=./advi --search-dir src:r=.. -- ./art_modern --qual_file_1 ../art/Illumina_profiles/HiSeq2500L125R1.txt --qual_file_2 ../art/Illumina_profiles/HiSeq2500L125R2.txt --i-file ../raw_data/ce11_chr1.fa --read_len 125 --mode wgs --lc se --i-parser memory --i-fcov 5 --parallel 0 --ins_rate_1 0.1 --del_rate_1 0.1 --o-sam test_small_se_wgs_memory.sam --pe_frag_dist_std_dev 20 --pe_frag_dist_mean 500

cd ..
