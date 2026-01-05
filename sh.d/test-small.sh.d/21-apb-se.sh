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
pigz -k "${OUT_DIR}"/out_se.fq
samtools view -h "${OUT_DIR}"/out_se.sam -o "${OUT_DIR}"/out_se.sam.bam
samtools view -h -u "${OUT_DIR}"/out_se.sam -o "${OUT_DIR}"/out_se.sam.uncompressed.bam
samtools view -T "${LAMBDA_PHAGE}" -h "${OUT_DIR}"/out_se.sam -o "${OUT_DIR}"/out_se.sam.cram
pigz -k "${OUT_DIR}"/out_se.sam
if [ -n "${WITH_NCBI_NGS:-}" ] && [ "${WITH_NCBI_NGS}" -eq 1 ]; then
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
    --i-file "${OUT_DIR}"/out_se.fq.gz \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_fq_gz.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
APB_EXEC \
    --i-format FASTQ \
    --i-file /dev/stdin \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_fq_cat.txt \
    --parallel "${PARALLEL}" \
    --old_behavior <"${OUT_DIR}"/out_se.fq
APB_EXEC \
    --i-format FASTQ \
    --i-file /dev/stdin \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_fq_gz_cat.txt \
    --parallel "${PARALLEL}" \
    --old_behavior <"${OUT_DIR}"/out_se.fq.gz
APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_sam.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
APB_EXEC \
    --i-format SAM \
    --i-file /dev/stdin \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_sam_cat.txt \
    --parallel "${PARALLEL}" \
    --old_behavior < "${OUT_DIR}"/out_se.sam
APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam.gz \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_sam_gz.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
APB_EXEC \
    --i-format SAM \
    --i-file /dev/stdin \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_sam_gz_cat.txt \
    --parallel "${PARALLEL}" \
    --old_behavior < "${OUT_DIR}"/out_se.sam.gz
APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam.bam \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_bam.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
APB_EXEC \
    --i-format BAM \
    --i-file /dev/stdin \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_bam_cat.txt \
    --parallel "${PARALLEL}" \
    --old_behavior < "${OUT_DIR}"/out_se.sam.bam
APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam.uncompressed.bam \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_ubam.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
APB_EXEC \
    --i-format BAM \
    --i-file /dev/stdin \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_ubam_cat.txt \
    --parallel "${PARALLEL}" \
    --old_behavior < "${OUT_DIR}"/out_se.sam.uncompressed.bam
APB_EXEC \
    --i-file "${OUT_DIR}"/out_se.sam.cram \
    --read_len "${RLEN}" \
    --o-file1 "${OUT_DIR}"/out_se_art_cxx_cram.txt \
    --parallel "${PARALLEL}" \
    --old_behavior
if [ -n "${WITH_NCBI_NGS:-}" ] && [ "${WITH_NCBI_NGS}" -eq 1 ]; then
    APB_EXEC \
        --i-file "${OUT_DIR}"/out_se.sra \
        --read_len "${RLEN}" \
        --o-file1 "${OUT_DIR}"/out_se_art_cxx_sra.txt \
        --parallel "${PARALLEL}" \
        --old_behavior
fi

cmp "${OUT_DIR}"/out_se_art_cxx_fq.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_fq_gz.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_fq_cat.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_fq_gz_cat.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_sam.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_sam_cat.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_sam_gz.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_sam_gz_cat.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_bam.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_bam_cat.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_ubam.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_ubam_cat.txt "${OUT_DIR}"/out_se_art_perl.txt
cmp "${OUT_DIR}"/out_se_art_cxx_cram.txt "${OUT_DIR}"/out_se_art_perl.txt

if [ -n "${WITH_NCBI_NGS:-}" ] && [ "${WITH_NCBI_NGS}" -eq 1 ]; then
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
