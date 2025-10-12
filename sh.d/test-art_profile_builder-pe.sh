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
    --o-sam "${TEST_DIR}"/out_pe.sam \
    --read_len 36 \
    --parallel 4 \
    --pe_frag_dist_mean 200 \
    --pe_frag_dist_std_dev 50 \
    --lc pe
samtools fastq \
    -1 "${TEST_DIR}"/out_pe.1.fq \
    -2 "${TEST_DIR}"/out_pe.2.fq \
    -N \
    "${TEST_DIR}"/out_pe.sam
# Generates: out_pe_art_perl_R1.txt out_pe_art_perl_R2.txt

deps/ART_profiler_illumina/art_profiler_illumina \
    "${TEST_DIR}"/out_pe_art_perl_ "${TEST_DIR}"/ fq

"${ART_MODERN_PATH}"/art_profile_builder \
    --i-file "${TEST_DIR}"/out_pe.1.fq \
    --read_len 36 \
    --o-file1 "${TEST_DIR}"/out_pe_art_cxx_fq_R1.txt \
    --parallel 2 \
    --i-num_threads 4 \
    --old_behavior
"${ART_MODERN_PATH}"/art_profile_builder \
    --i-file "${TEST_DIR}"/out_pe.2.fq \
    --read_len 36 \
    --o-file1 "${TEST_DIR}"/out_pe_art_cxx_fq_R2.txt \
    --parallel 2 \
    --i-num_threads 4 \
    --old_behavior
"${ART_MODERN_PATH}"/art_profile_builder \
    --i-file "${TEST_DIR}"/out_pe.sam \
    --read_len 36 \
    --is_pe \
    --o-file1 "${TEST_DIR}"/out_pe_art_cxx_sam_R1.txt \
    --o-file2 "${TEST_DIR}"/out_pe_art_cxx_sam_R2.txt \
    --parallel 2 \
    --i-num_threads 4 \
    --old_behavior

cmp "${TEST_DIR}"/out_pe_art_cxx_fq_R1.txt "${TEST_DIR}"/out_pe_art_perl_R1.txt
cmp "${TEST_DIR}"/out_pe_art_cxx_fq_R2.txt "${TEST_DIR}"/out_pe_art_perl_R2.txt

cmp "${TEST_DIR}"/out_pe_art_cxx_sam_R1.txt "${TEST_DIR}"/out_pe_art_perl_R1.txt
cmp "${TEST_DIR}"/out_pe_art_cxx_sam_R2.txt "${TEST_DIR}"/out_pe_art_perl_R2.txt
