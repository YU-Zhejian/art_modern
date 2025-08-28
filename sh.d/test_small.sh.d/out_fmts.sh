# shellcheck shell=bash

"${ART}" \
    --qual_file_1 data/Illumina_profiles/HiSeq2500L150R1.txt \
    --sep_flag \
    --i-file "${MRNA_HEAD}" \
    --read_len 150 \
    --mode template \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq "${OUT_DIR}"/test_small_se_template_memory_sep.fastq \
    --o-fasta "${OUT_DIR}"/test_small_se_template_memory_sep.fasta \
    --o-pwa "${OUT_DIR}"/test_small_se_template_memory_sep.pwa \
    --o-sam "${OUT_DIR}"/test_small_se_template_memory_sep.sam \
    --o-hl_sam "${OUT_DIR}"/test_small_se_template_memory_sep.hl.bam \
    --o-hl_sam-num_threads 2 \
    --o-hl_sam-compress_level u \
    --o-hl_sam-write_bam
if type -p fastqc &> /dev/null && type -p x-www-browser &>/dev/null; then
    fastqc "${OUT_DIR}"/test_small_se_template_memory_sep.fastq
    # Open the browser and ignore what's happening afterwards
    {
        x-www-browser "${OUT_DIR}"/test_small_se_template_memory_sep_fastqc.html &>/dev/null || true
    } &
    sleep 3
fi
if [ "${FORMAT_ONLY:-}" = "1" ]; then
    exit 0
fi
rm -fr "${OUT_DIR}"/test_small_se_template_memory_sep.* "${OUT_DIR}"/test_small_se_template_memory_sep_fastqc.*
