# shellcheck shell=bash

AM_EXEC \
    --i-file "${CE11_CHR1}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --i-seed 42 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq "${OUT_DIR}"/test_small_se_wgs_memory_sep_1.fastq
merge_file "${OUT_DIR}"/test_small_se_wgs_memory_sep_1.fastq
cat "${OUT_DIR}"/test_small_se_wgs_memory_sep_1.fastq |
    seqkit sort >"${OUT_DIR}"/test_small_se_wgs_memory_sep_1.fastq.srt.fq

AM_EXEC \
    --i-file "${CE11_CHR1}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --i-seed 42 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq
merge_file "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq
cat "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq |
    seqkit sort >"${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq.srt.fq
cmp \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep_1.fastq.srt.fq \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq.srt.fq

AM_EXEC \
    --i-file "${CE11_CHR1}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --i-seed 43 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq
merge_file "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq
cat "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq |
    seqkit sort >"${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq.srt.fq

! cmp \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep_1.fastq.srt.fq \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.fastq.srt.fq &>/dev/null

rm -fr \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep_1.* \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep_2.*

assert_cleandir
