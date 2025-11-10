# shellcheck shell=bash

FCOV=0.2
for parser in memory htslib; do
    for lc in se pe mp; do
        AM_EXEC \
            --i-file "${CE11_CHR1}" \
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
        python sh.d/test-small.sh.d/validate_rlen.py \
            "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
            150
        rm -fr "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
        assert_cleandir
    done
done

# We now test whether equal --read_len_1 and --read_len_2 is working.
parser=memory
for lc in pe mp; do
    AM_EXEC \
        --i-file "${CE11_CHR1}" \
        --i-batch_size 100 \
        --mode wgs \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --read_len 100 \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 500
    merge_file "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
    python sh.d/test-small.sh.d/validate_rlen.py \
        "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
        100
    rm -fr "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
    assert_cleandir

    AM_EXEC \
        --i-file "${CE11_CHR1}" \
        --i-batch_size 100 \
        --mode wgs \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --read_len_1 100 \
        --read_len_2 100 \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 500
    merge_file "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
    python sh.d/test-small.sh.d/validate_rlen.py \
        "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
        100
    rm -fr "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
    assert_cleandir
done

lc=se
AM_EXEC \
    --i-file "${CE11_CHR1}" \
    --i-batch_size 100 \
    --mode wgs \
    --lc "${lc}" \
    --i-parser "${parser}" \
    --i-fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --read_len 100 \
    --o-fastq "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500
merge_file "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
python sh.d/test-small.sh.d/validate_rlen.py \
    "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
    100
rm -fr "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
assert_cleandir

AM_EXEC \
    --i-file "${CE11_CHR1}" \
    --i-batch_size 100 \
    --mode wgs \
    --lc "${lc}" \
    --i-parser "${parser}" \
    --i-fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --read_len_1 100 \
    --o-fastq "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500
merge_file "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
python sh.d/test-small.sh.d/validate_rlen.py \
    "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq \
    100
rm -fr "${OUT_DIR}"/test_small_"${lc}"_wgs_"${parser}".fq
assert_cleandir

unset FCOV
