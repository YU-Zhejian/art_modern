# shellcheck shell=bash

FCOV=10
parser=memory
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file "${MRNA_HEAD}" \
        --read_len 125 \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".sam \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}" "${MRNA_HEAD}"
    python sh.d/test-small.sh.d/validate_cov.py \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq \
        "${MRNA_HEAD}" \
        "${FCOV}" \
        CONST_COV \
        IS_TEMPLATE
done

parser=stream
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file "${MRNA_HEAD}" \
        --read_len 125 \
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
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".hl "${MRNA_HEAD}"
    python sh.d/test-small.sh.d/validate_cov.py \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".fq \
        "${MRNA_HEAD}" \
        "${FCOV}" \
        CONST_COV \
        IS_TEMPLATE
done
rm -fr "${OUT_DIR}"/test_small_??_template_stream.hl.sam
unset FCOV
