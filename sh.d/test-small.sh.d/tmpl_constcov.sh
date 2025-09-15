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
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".sam
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}" "${MRNA_HEAD}"
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
        --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".hl.sam
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}".hl "${MRNA_HEAD}"
done
rm -fr "${OUT_DIR}"/test_small_??_template_stream.hl.sam
unset FCOV
