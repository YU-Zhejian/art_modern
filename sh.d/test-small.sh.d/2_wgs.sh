# shellcheck shell=bash

FCOV=0.2
for parser in memory htslib; do
    for lc in se pe mp; do
        AM_EXEC \
            --builtin_qual_file HiSeq2500_125bp \
            --i-file "${CE11_CHR1}" \
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
        merge_file "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
        merge_file "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".sam
        sam2bam "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}" "${CE11_CHR1}"
        python sh.d/test-small.sh.d/validate_cov.py \
            "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
            "${CE11_CHR1}" \
            "${FCOV}" \
            CONST_COV \
            NOT_TEMPLATE
        rm -fr "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
        assert_cleandir
    done
done
unset FCOV
