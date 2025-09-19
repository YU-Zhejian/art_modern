# shellcheck shell=bash

for coverage in stranded strandless; do
    parser=memory
    for lc in se pe mp; do
        "${ART}" \
            --builtin_qual_file HiSeq2500_125bp \
            --i-file "${MRNA_HEAD}" \
            --read_len 125 \
            --mode template \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov data/raw_data/ce11.mRNA_head.cov_"${coverage}".tsv \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".sam
        sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}" "${MRNA_HEAD}"
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
            --i-fcov data/raw_data/ce11.mRNA_head.cov_"${coverage}".tsv \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl.sam
        sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl "${MRNA_HEAD}"
    done
    rm -fr "${OUT_DIR}"/test_small_??_template_"${parser}"_"${coverage}".hl.sam
done
