# shellcheck shell=bash
# Test SAM Tags

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-sam "${OUT_DIR}"/test_small_se_wgs_memory.sam \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam \
    --o-hl_sam-num_threads 2 \
    --o-hl_sam-compress_level u \
    --o-hl_sam-write_bam

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.sam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.sam NM MD
validate_read_quals "${OUT_DIR}"/test_small_se_wgs_memory.sam true
sam2bam "${OUT_DIR}"/test_small_se_wgs_memory "${LAMBDA_PHAGE}"

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam OA NM MD
validate_read_quals "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam true

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-sam "${OUT_DIR}"/test_small_se_wgs_memory.sam \
    --o-sam-without_tag_NM

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.sam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.sam MD
sam2bam "${OUT_DIR}"/test_small_se_wgs_memory "${LAMBDA_PHAGE}"

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-sam "${OUT_DIR}"/test_small_se_wgs_memory.sam \
    --o-sam-without_tag_MD

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.sam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.sam NM
sam2bam "${OUT_DIR}"/test_small_se_wgs_memory "${LAMBDA_PHAGE}"

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-sam "${OUT_DIR}"/test_small_se_wgs_memory.sam \
    --o-sam-without_tag_MD \
    --o-sam-without_tag_NM

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.sam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.sam
sam2bam "${OUT_DIR}"/test_small_se_wgs_memory "${LAMBDA_PHAGE}"

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam \
    --o-hl_sam-without_tag_NM \
    --o-hl_sam-write_bam

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam MD OA

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam \
    --o-hl_sam-without_tag_MD \
    --o-hl_sam-write_bam

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam NM OA

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam \
    --o-hl_sam-without_tag_MD \
    --o-hl_sam-without_tag_NM \
    --o-hl_sam-write_bam

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam OA

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam \
    --o-hl_sam-without_tag_NM \
    --o-hl_sam-without_tag_OA \
    --o-hl_sam-write_bam

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam MD

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam \
    --o-hl_sam-without_tag_MD \
    --o-hl_sam-without_tag_OA \
    --o-hl_sam-write_bam

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam NM

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam \
    --o-hl_sam-without_tag_MD \
    --o-hl_sam-without_tag_NM \
    --o-hl_sam-without_tag_OA \
    --o-hl_sam-write_bam

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam
validate_read_tags "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam

AM_EXEC \
    --i-file "${LAMBDA_PHAGE}" \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --o-sam "${OUT_DIR}"/test_small_se_wgs_memory.sam \
    --o-sam-no_qual \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam \
    --o-hl_sam-no_qual \
    --o-hl_sam-write_bam

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.sam
validate_read_quals "${OUT_DIR}"/test_small_se_wgs_memory.sam false

merge_file "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam
validate_read_quals "${OUT_DIR}"/test_small_se_wgs_memory.hl.bam false

rm -fr "${OUT_DIR}"/test_small_se_wgs_memory.*
assert_cleandir
