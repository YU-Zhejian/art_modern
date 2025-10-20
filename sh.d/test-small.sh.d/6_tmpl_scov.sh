# shellcheck shell=bash

for coverage in stranded strandless; do
    parser=memory
    for lc in se pe mp; do
        "${ART_CMD_ASSEMBLED[@]}" \
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
            --o-sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".sam \
            --o-fastq "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq
        merge_file "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".sam
        merge_file "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq
        sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}" "${MRNA_HEAD}"
        python sh.d/test-small.sh.d/validate_cov.py \
            "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq \
            "${MRNA_HEAD}" \
            data/raw_data/ce11.mRNA_head.cov_"${coverage}".tsv \
            COV_TSV \
            IS_TEMPLATE
        rm -fr "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq
        assert_cleandir
    done

    parser=stream
    for lc in se pe mp; do
        "${ART_CMD_ASSEMBLED[@]}" \
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
            --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl.sam \
            --o-fastq "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq
        merge_file "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl.sam
        merge_file "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq
        sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl "${MRNA_HEAD}"
        python sh.d/test-small.sh.d/validate_cov.py \
            "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq \
            "${MRNA_HEAD}" \
            data/raw_data/ce11.mRNA_head.cov_"${coverage}".tsv \
            COV_TSV \
            IS_TEMPLATE
        rm -fr \
            "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq \
            "${OUT_DIR}"/test_small_??_template_"${parser}"_"${coverage}".hl.sam
        assert_cleandir
    done
done
