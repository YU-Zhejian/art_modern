# shellcheck shell=bash

FCOV=10
parser=memory
for lc in se pe mp; do
    AM_EXEC \
        --i-file "${MRNA_HEAD}" \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".sam \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq
    merge_file "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".sam
    merge_file "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq
    validate_template \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".sam "${lc}"
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}" "${MRNA_HEAD}"
    validate_cov \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq \
        "${MRNA_HEAD}" \
        "${FCOV}" \
        CONST_COV \
        IS_TEMPLATE
    rm -fr "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq
    assert_cleandir
done

parser=stream
for lc in se pe mp; do
    AM_EXEC \
        --i-file "${MRNA_HEAD}" \
        --i-batch_size 100 \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".hl.sam \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq
    merge_file "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".hl.sam
    merge_file "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".hl "${MRNA_HEAD}"
    validate_cov \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq \
        "${MRNA_HEAD}" \
        "${FCOV}" \
        CONST_COV \
        IS_TEMPLATE
    rm -fr \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".hl.sam
    assert_cleandir
done
unset FCOV
