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
merge_file "${OUT_DIR}"/out_se.sam
samtools fastq \
    -0 "${OUT_DIR}"/out_se.fq \
    -n \
    "${OUT_DIR}"/out_se.sam
if [ -n ${WITH_SRA:-} ] && [ "${WITH_SRA}" -eq 1 ]; then
    latf-load \
        --no-readnames \
        -p ILLUMINA \
        -o "${OUT_DIR}"/out_se.sra.d \
        -L info \
        -q PHRED_33 \
        "${OUT_DIR}"/out_se.fq
    kar -f -c "${OUT_DIR}"/out_se.sra -d "${OUT_DIR}"/out_se.sra.d --md5
fi

art_profile_illumina "${OUT_DIR}"/out_se_art_perl "${OUT_DIR}"/ fq

APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.fq \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_fq.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_sam.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
if [ -n ${WITH_SRA:-} ] && [ "${WITH_SRA}" -eq 1 ]; then
    APB_EXEC \
        --i-file "${OUT_DIR}"/out_se.sra \
        --read_len "${RLEN}" \
        --o-file1 "${OUT_DIR}"/out_se_art_cxx_sra.txt \
        --parallel "${PARALLEL}" \
        --old_behavior
fi

cmp "${OUT_DIR}"/out_se_art_cxx_fq.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_sam.txt "${OUT_DIR}"/out_se_art_perl.txt
if [ -n ${WITH_SRA:-} ] && [ "${WITH_SRA}" -eq 1 ]; then
    cmp "${OUT_DIR}"/out_se_art_cxx_sra.txt "${OUT_DIR}"/out_se_art_perl.txt
fi

APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam \
    --read_len_1 "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_sam.txt \
    --parallel "${PARALLEL}" \
    --old_behavior

if [ ! "${SET_X:-0}" == "1" ]; then
    set -x
fi
[ "$(grep -c -v '^$' <"${OUT_DIR}"/out_se_art_cxx_sam.txt)" == "432" ]
if [ ! "${SET_X:-0}" == "1" ]; then
    set +x
fi

APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam \
    --read_len_1 10 \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_sam.txt \
    --parallel "${PARALLEL}" \
    --old_behavior

if [ ! "${SET_X:-0}" == "1" ]; then
    set -x
fi
[ "$(grep -c -v '^$' <"${OUT_DIR}"/out_se_art_cxx_sam.txt)" == "120" ]
if [ ! "${SET_X:-0}" == "1" ]; then
    set +x
fi

rm -fr "${OUT_DIR:?}"/*
assert_cleandir
unset RLEN
