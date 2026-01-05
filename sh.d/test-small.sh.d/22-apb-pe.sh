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
merge_file "${OUT_DIR}"/out_pe.sam
samtools fastq \
    -1 "${OUT_DIR}"/out_pe.1.fq \
    -2 "${OUT_DIR}"/out_pe.2.fq \
    -N \
    "${OUT_DIR}"/out_pe.sam

if [ -n "${WITH_NCBI_NGS:-}" ] && [ "${WITH_NCBI_NGS}" -eq 1 ]; then
    latf-load \
        --no-readnames \
        -p ILLUMINA \
        -o "${OUT_DIR}"/out_pe.sra.d \
        -L info \
        -q PHRED_33 \
        "${OUT_DIR}"/out_pe.1.fq "${OUT_DIR}"/out_pe.2.fq
    kar -f -c "${OUT_DIR}"/out_pe.sra -d "${OUT_DIR}"/out_pe.sra.d --md5
fi

art_profile_illumina "${OUT_DIR}"/out_pe_art_perl_ "${OUT_DIR}"/ fq

APB_EXEC \
    --i-file "${OUT_DIR}"/out_pe.1.fq \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_pe_art_cxx_fq_R1.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
APB_EXEC \
    --i-file "${OUT_DIR}"/out_pe.2.fq \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_pe_art_cxx_fq_R2.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
APB_EXEC \
    --i-file "${OUT_DIR}"/out_pe.sam \
    --read_len "${RLEN}" \
    --is_pe \
    --o-file1 "${OUT_DIR}"/out_pe_art_cxx_sam_R1.txt \
    --o-file2 "${OUT_DIR}"/out_pe_art_cxx_sam_R2.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
if [ -n "${WITH_NCBI_NGS:-}" ] && [ "${WITH_NCBI_NGS}" -eq 1 ]; then
    APB_EXEC \
        --i-file "${OUT_DIR}"/out_pe.sra \
        --read_len "${RLEN}" \
        --is_pe \
        --o-file1 "${OUT_DIR}"/out_pe_art_cxx_sra_R1.txt \
        --o-file2 "${OUT_DIR}"/out_pe_art_cxx_sra_R2.txt \
        --parallel "${PARALLEL}" \
        --old_behavior
fi

cmp "${OUT_DIR}"/out_pe_art_cxx_fq_R1.txt "${OUT_DIR}"/out_pe_art_perl_R1.txt
cmp "${OUT_DIR}"/out_pe_art_cxx_fq_R2.txt "${OUT_DIR}"/out_pe_art_perl_R2.txt

cmp "${OUT_DIR}"/out_pe_art_cxx_sam_R1.txt "${OUT_DIR}"/out_pe_art_perl_R1.txt
cmp "${OUT_DIR}"/out_pe_art_cxx_sam_R2.txt "${OUT_DIR}"/out_pe_art_perl_R2.txt

if [ -n "${WITH_NCBI_NGS:-}" ] && [ "${WITH_NCBI_NGS}" -eq 1 ]; then
    cmp "${OUT_DIR}"/out_pe_art_cxx_sra_R1.txt "${OUT_DIR}"/out_pe_art_perl_R1.txt
    cmp "${OUT_DIR}"/out_pe_art_cxx_sra_R2.txt "${OUT_DIR}"/out_pe_art_perl_R2.txt
fi

APB_EXEC \
    --i-file "${OUT_DIR}"/out_pe.sam \
    --read_len_1 "${RLEN}" \
    --read_len_2 "${RLEN}" \
    --is_pe \
    --o-file1 "${OUT_DIR}"/out_pe_art_cxx_sam_R1.txt \
    --o-file2 "${OUT_DIR}"/out_pe_art_cxx_sam_R2.txt \
    --parallel "${PARALLEL}" \
    --i-num_threads 4 \
    --old_behavior

if [ ! "${SET_X:-0}" == "1" ]; then
    set -x
fi
[ "$(grep -c -v '^$' <"${OUT_DIR}"/out_pe_art_cxx_sam_R1.txt)" == "432" ]
[ "$(grep -c -v '^$' <"${OUT_DIR}"/out_pe_art_cxx_sam_R2.txt)" == "432" ]
if [ ! "${SET_X:-0}" == "1" ]; then
    set +x
fi

APB_EXEC \
    --i-file "${OUT_DIR}"/out_pe.sam \
    --read_len_1 10 \
    --read_len_2 "${RLEN}" \
    --is_pe \
    --o-file1 "${OUT_DIR}"/out_pe_art_cxx_sam_R1.txt \
    --o-file2 "${OUT_DIR}"/out_pe_art_cxx_sam_R2.txt \
    --parallel "${PARALLEL}" \
    --i-num_threads 4 \
    --old_behavior

if [ ! "${SET_X:-0}" == "1" ]; then
    set -x
fi
[ "$(grep -c -v '^$' <"${OUT_DIR}"/out_pe_art_cxx_sam_R1.txt)" == "120" ]
[ "$(grep -c -v '^$' <"${OUT_DIR}"/out_pe_art_cxx_sam_R2.txt)" == "432" ]
if [ ! "${SET_X:-0}" == "1" ]; then
    set +x
fi

rm -fr "${OUT_DIR:?}"/*
assert_cleandir

unset RLEN
