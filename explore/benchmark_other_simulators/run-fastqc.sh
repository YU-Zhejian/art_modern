#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
# shellcheck disable=SC2155
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

export OUT_DIR="$(pwd)"/data_out
export ART_MODERN_THREADS=20

mkdir -p "${OUT_DIR}"

bin/wgsim \
    -1 150 -2 150 -N 1022922 -d 300 -s 20 -r 0 \
    data/yeast.fa \
    "${OUT_DIR}"/yeast_wgsim_1.fq "${OUT_DIR}"/yeast_wgsim_2.fq &

bin/dwgsim \
    -1 150 -2 150 -C 3 -d 300 -s 20 -o 1 -r 0 -y 0 \
    data/yeast.fa \
    "${OUT_DIR}"/yeast_dwgsim &

bin/art_original \
    --in data/yeast.fa --out "${OUT_DIR}"/yeast_art_ \
    -f 3 --len 150 --mflen 300 --sdev 20 --noALN --paired --seqSys HS25 &

wait
opt/art_modern_build/art_modern \
    --mode wgs --lc pe \
    --i-file data/yeast.fa --i-fcov 3 --read_len 150 \
    --o-fastq "${OUT_DIR}"/yeast_art_modern_wgs_memory.fastq \
    --qual_file_1 ../../data/Illumina_profiles/HiSeq2500L150R1filter.txt \
    --qual_file_2 ../../data/Illumina_profiles/HiSeq2500L150R2filter.txt \
    --pe_frag_dist_mean 300 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"
seqtk seq -1 "${OUT_DIR}"/yeast_art_modern_wgs_memory.fastq >"${OUT_DIR}"/yeast_art_modern_wgs_memory_1.fastq
seqtk seq -2 "${OUT_DIR}"/yeast_art_modern_wgs_memory.fastq >"${OUT_DIR}"/yeast_art_modern_wgs_memory_2.fastq
rm -f "${OUT_DIR}"/yeast_art_modern_wgs_memory.fastq

fastqc \
    --extract \
    --svg \
    --threads 20 \
    "${OUT_DIR}"/*.fastq \
    "${OUT_DIR}"/*.fq
rm -fr "${OUT_DIR}"/*.fastq "${OUT_DIR}"/*.fq
