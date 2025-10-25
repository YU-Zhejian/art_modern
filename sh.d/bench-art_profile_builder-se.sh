#!/usr/bin/env bash
set -uex
ART_MODERN_PATH="opt/build_release_install/bin"
TEST_DIR="opt/test-art_profile_builder"
rm -rf "${TEST_DIR}"
mkdir -p "${TEST_DIR}"

"${ART_MODERN_PATH}"/art_modern \
    --mode wgs \
    --i-file data/raw_data/lambda_phage.fa \
    --i-fcov 400000 \
    --builtin_qual_file HiSeq1000_100bp \
    --o-fastq "${TEST_DIR}"/out_se.fq \
    --read_len 100 \
    --parallel 16 \
    --lc se
pigz -k -9 -f -p16 opt/test-art_profile_builder/out_se.fq -vv

# real    53m17.432s
# user    53m0.175s
# sys     0m16.992s
time deps/ART_profiler_illumina/art_profiler_illumina \
    "${TEST_DIR}"/out_se_art_perl_uncompressed "${TEST_DIR}"/ fq 16

time deps/ART_profiler_illumina/art_profiler_illumina \
    "${TEST_DIR}"/out_se_art_perl_compressed "${TEST_DIR}"/ fq.gz 16

# SAM/BAM parser:
#     real    1m28.444s
#     user    6m25.348s
#     sys     1m31.002s
# FASTQ parser:
#     real    1m27.234s
#     user    3m42.257s
#     sys     1m46.697s
time "${ART_MODERN_PATH}"/art_profile_builder \
    --i-file "${TEST_DIR}"/out_se.fq \
    --read_len 36 \
    --o-file1 "${TEST_DIR}"/out_se_art_cxx_new_uncompressed_fq.txt \
    --parallel 8 \
    --i-num_threads 4 \
    --old_behavior

# SAM/BAM parser:
#     real    4m7.237s
#     user    31m8.373s
#     sys     0m24.162s
# FASTQ parser:
#     real    6m25.944s
#     user    46m18.036s
#     sys     0m33.816s
time "${ART_MODERN_PATH}"/art_profile_builder \
    --i-file "${TEST_DIR}"/out_se.fq.gz \
    --read_len 36 \
    --o-file1 "${TEST_DIR}"/out_se_art_cxx_new_compressed_fq.txt \
    --parallel 8 \
    --i-num_threads 4 \
    --old_behavior

# Rebuilt using HTSLib parser
