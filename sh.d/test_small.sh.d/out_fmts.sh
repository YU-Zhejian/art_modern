# shell=bash

"${ART}" \
    --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt \
    --sep_flag \
    --i-file "${MRNA_HEAD}" \
    --read_len 125 \
    --mode template \
    --lc se \
    --i-parser memory \
    --i-fcov 1 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq "${OUT_DIR}"/test_small_se_template_memory_sep.fastq \
    --o-pwa "${OUT_DIR}"/test_small_se_template_memory_sep.pwa \
    --o-sam "${OUT_DIR}"/test_small_se_template_memory_sep.sam \
    --o-hl_sam "${OUT_DIR}"/test_small_se_template_memory_sep.hl.bam \
    --o-hl_sam-write_bam
rm -fr "${OUT_DIR}"/test_small_se_template_memory_sep.*
