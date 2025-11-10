# shellcheck shell=bash

RLEN=36
AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --i-fcov 0.4 \
    --builtin_qual_file GA1_36bp \
    --o-sam "${OUT_DIR}"/out_se.sam \
    --read_len "${RLEN}" \
    --parallel "${PARALLEL}" \
    --lc se
samtools fastq \
    -0 "${OUT_DIR}"/out_se.fq \
    -n \
    "${OUT_DIR}"/out_se.sam

deps/ART_profiler_illumina/art_profiler_illumina \
    "${OUT_DIR}"/out_se_art_perl "${OUT_DIR}"/ fq

ABP_EXEC \
    --i-file "${OUT_DIR}"/out_se.fq \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_fq.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
ABP_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_sam.txt \
    --parallel "${PARALLEL}" \
    --old_behavior

cmp "${OUT_DIR}"/out_se_art_cxx_fq.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_sam.txt "${OUT_DIR}"/out_se_art_perl.txt

rm -fr "${OUT_DIR}"/*
assert_cleandir
unset RLEN
