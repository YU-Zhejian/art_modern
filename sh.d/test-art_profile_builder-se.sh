#!/usr/bin/env bash
set -uex
ART_MODERN_PATH="opt/build_release_install/bin"
TEST_DIR="opt/test-art_profile_builder"
rm -rf "${TEST_DIR}"
mkdir -p "${TEST_DIR}"

"${ART_MODERN_PATH}"/art_modern \
    --mode wgs \
    --i-file data/raw_data/lambda_phage.fa \
    --i-fcov 0.4 \
    --builtin_qual_file GA1_36bp \
    --o-sam "${TEST_DIR}"/out_se.sam \
    --read_len 36 \
    --parallel 4 \
    --lc se
samtools fastq \
    -0 "${TEST_DIR}"/out_se.fq \
    -n \
    "${TEST_DIR}"/out_se.sam
# Generates: out_se_art_perl_R1.txt

deps/ART_profiler_illumina/art_profiler_illumina \
    "${TEST_DIR}"/out_se_art_perl "${TEST_DIR}"/ fq

"${ART_MODERN_PATH}"/art_profile_builder \
    --i-file "${TEST_DIR}"/out_se.fq \
    --read_len 36 \
    --o-file1 "${TEST_DIR}"/out_se_art_cxx_fq.txt \
    --parallel 2 \
    --i-num_threads 4 \
    --old_behavior
"${ART_MODERN_PATH}"/art_profile_builder \
    --i-file "${TEST_DIR}"/out_se.sam \
    --read_len 36 \
    --o-file1 "${TEST_DIR}"/out_se_art_cxx_sam.txt \
    --parallel 2 \
    --i-num_threads 4 \
    --old_behavior

cmp "${TEST_DIR}"/out_se_art_cxx_fq.txt "${TEST_DIR}"/out_se_art_perl.txt
cmp "${TEST_DIR}"/out_se_art_cxx_sam.txt "${TEST_DIR}"/out_se_art_perl.txt
