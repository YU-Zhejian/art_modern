# shellcheck shell=bash

RLEN=36
AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --i-fcov 0.4 \
    --builtin_qual_file GA1_36bp \
    --o-sam "${OUT_DIR}"/out_pe.sam \
    --read_len "${RLEN}" \
    --parallel "${PARALLEL}" \
    --pe_frag_dist_mean 200 \
    --pe_frag_dist_std_dev 50 \
    --lc pe
samtools fastq \
    -1 "${OUT_DIR}"/out_pe.1.fq \
    -2 "${OUT_DIR}"/out_pe.2.fq \
    -N \
    "${OUT_DIR}"/out_pe.sam

deps/ART_profiler_illumina/art_profiler_illumina \
    "${OUT_DIR}"/out_pe_art_perl_ "${OUT_DIR}"/ fq

ABP_EXEC \
    --i-file "${OUT_DIR}"/out_pe.1.fq \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_pe_art_cxx_fq_R1.txt \
    --parallel "${PARALLEL}" \
    --i-num_threads 4 \
    --old_behavior
ABP_EXEC \
    --i-file "${OUT_DIR}"/out_pe.2.fq \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_pe_art_cxx_fq_R2.txt \
    --parallel "${PARALLEL}" \
    --i-num_threads 4 \
    --old_behavior
ABP_EXEC \
    --i-file "${OUT_DIR}"/out_pe.sam \
    --read_len "${RLEN}" \
    --is_pe \
    --o-file1 "${OUT_DIR}"/out_pe_art_cxx_sam_R1.txt \
    --o-file2 "${OUT_DIR}"/out_pe_art_cxx_sam_R2.txt \
    --parallel "${PARALLEL}" \
    --i-num_threads 4 \
    --old_behavior

cmp "${OUT_DIR}"/out_pe_art_cxx_fq_R1.txt "${OUT_DIR}"/out_pe_art_perl_R1.txt
cmp "${OUT_DIR}"/out_pe_art_cxx_fq_R2.txt "${OUT_DIR}"/out_pe_art_perl_R2.txt

cmp "${OUT_DIR}"/out_pe_art_cxx_sam_R1.txt "${OUT_DIR}"/out_pe_art_perl_R1.txt
cmp "${OUT_DIR}"/out_pe_art_cxx_sam_R2.txt "${OUT_DIR}"/out_pe_art_perl_R2.txt

rm -fr "${OUT_DIR}"/*
assert_cleandir
