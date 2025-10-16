# shellcheck shell=bash

FCOV=0.2
for parser in memory htslib; do
    for lc in se pe mp; do
        "${ART}" \
            --builtin_qual_file HiSeq2500_125bp \
            --i-file data/raw_data/ce11_chr1.fa \
            --read_len 125 \
            --i-batch_size 100 \
            --mode wgs \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov "${FCOV}" \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-sam "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".sam \
            --o-fastq "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
            --pe_frag_dist_std_dev 20 \
            --pe_frag_dist_mean 500
        sam2bam "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}" data/raw_data/ce11_chr1.fa
        python sh.d/test-small.sh.d/validate_cov.py \
            "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
            data/raw_data/ce11_chr1.fa \
            "${FCOV}" \
            CONST_COV \
            NOT_TEMPLATE
    done
done
unset FCOV
