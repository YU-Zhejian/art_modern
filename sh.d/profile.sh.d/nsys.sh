#!/usr/bin/env bash
set -ue
export NVIDIA_HPC_SDK_PATH="${NVIDIA_HPC_SDK_PATH:-/opt/nvidia/hpc_sdk/Linux_x86_64/2024/}"

if [ ! -d "${NVIDIA_HPC_SDK_PATH}" ]; then
    echo "NVIDIA HPC SDK not found at ${NVIDIA_HPC_SDK_PATH}"
    exit 1
else
    echo "NVIDIA HPC SDK found at ${NVIDIA_HPC_SDK_PATH}"
fi

PROFILE_DIR=opt/build_profile_nsys
mkdir -p "${PROFILE_DIR}"

env -C "${PROFILE_DIR}" PATH="${NVIDIA_HPC_SDK_PATH}/compilers/bin/:${PATH:-}" cmake \
    -DCMAKE_C_COMPILER=nvc \
    -DCMAKE_CXX_COMPILER=nvc++ \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DUSE_RANDOM_GENERATOR=GSL \
    -G Ninja "$(pwd)"

env -C "${PROFILE_DIR}" PATH="${NVIDIA_HPC_SDK_PATH}/compilers/bin/:${PATH:-}" ninja -j120

"${PROFILE_DIR}"/art_modern --version
ldd "${PROFILE_DIR}"/art_modern

sudo nsys profile \
    --trace=openmp \
    --force-overwrite=true \
    --output="${PROFILE_DIR}"/report1.nsys-rep \
    ./art_modern \
    --builtin_qual_file HiSeq2500_125bp \
    --i-file data/raw_data/ce11_chr1.fa \
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
    --o-sam "${PROFILE_DIR}"/test_small_se_wgs_memory.bam \
    --o-sam-write_bam \
    --o-hl_sam "${PROFILE_DIR}"/test_small_se_wgs_memory.hl.sam \
    --o-fastq "${PROFILE_DIR}"/test_small_se_wgs_memory.fastq \
    --o-pwa "${PROFILE_DIR}"/test_small_se_wgs_memory.pwa
nsys-ui "${PROFILE_DIR}"/report1.nsys-rep
cd ..
