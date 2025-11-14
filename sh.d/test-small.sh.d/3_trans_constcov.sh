# shellcheck shell=bash

# Trans mode with constant coverage
parser=memory
FCOV=10
for lc in se pe mp; do
    AM_EXEC \
        --i-file "${MRNA_HEAD}" \
        --mode trans \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}".sam \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}".fq \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 500
    merge_file "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}".sam
    merge_file "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}".fq
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}" "${MRNA_HEAD}"
    validate_cov \
        "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}".fq \
        "${MRNA_HEAD}" \
        "${FCOV}" \
        CONST_COV \
        NOT_TEMPLATE
    rm -fr "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}".fq
    assert_cleandir
done
# No need to test stream FASTA parser
unset FCOV
