#!/usr/bin/env bash
export NVIDIA_HPC_SDK_PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/2024/"

set -ue
rm -fr build_profile
mkdir -p build_profile
cd build_profile
env PATH="${NVIDIA_HPC_SDK_PATH}/compilers/bin/:${PATH:-}" cmake \
    -DCMAKE_C_COMPILER=nvc \
    -DCMAKE_CXX_COMPILER=nvc++ \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=GSL \
    -DUSE_HTSLIB=hts \
    -G Ninja .. \
    -DCMAKE_PREFIX_PATH=/home/yuzj/opt/samtools-1.21/ \
    -DCMAKE_INCLUDE_PATH=/home/yuzj/opt/samtools-1.21/include

ninja -j120

./art_modern --version
ldd ./art_modern

sudo nsys profile \
    --trace=openmp \
    --force-overwrite=true \
    --output=report1.nsys-rep \
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
/usr/lib/nsight-systems/host-linux-x64/QdstrmImporter -i report1.qdstrm -o report1.nsys-rep
nsys-ui report1.nsys-rep
cd ..
#
#    --o-sam test_small_se_wgs_memory.sam \
#    --o-sam-write_bam \
