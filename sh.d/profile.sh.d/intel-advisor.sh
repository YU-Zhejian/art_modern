#!/usr/bin/env bash
# shellcheck disable=SC2317
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
rm -fr build_profile
mkdir -p build_profile
cd build_profile
cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=ONEMKL \
    -DUSE_HTSLIB=hts \
    -G Ninja .. \
    -DCMAKE_PREFIX_PATH=/home/yuzj/opt/samtools-1.21/ \
    -DCMAKE_INCLUDE_PATH=/home/yuzj/opt/samtools-1.21/include

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
    --lc pe \
    --i-parser memory \
    --i-fcov 4 \
    --parallel 0 \
    --ins_rate_1 0.1 \
    --del_rate_1 0.1 \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500 \
    --o-sam test_small_se_wgs_memory.bam \
    --o-sam-write_bam \
    --o-hl_sam test_small_se_wgs_memory.hl.sam \
    --o-fastq test_small_se_wgs_memory.fastq \
    --o-pwa test_small_se_wgs_memory.pwa
advisor-gui ./advi_results
cd ..
#
#    --o-sam test_small_se_wgs_memory.sam \
#    --o-sam-write_bam \
