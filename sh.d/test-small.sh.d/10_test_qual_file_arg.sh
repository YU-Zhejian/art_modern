#shellcheck shell=bash
export KEEPLOG=1
touch "${OUT_DIR}"/devnull

# This test is added to ensure the quality profile args are working
AM_EXEC \
    --i-file "${OUT_DIR}"/devnull \
    --lc "pe" \
    --i-parser memory \
    --i-fcov 0 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 200 \
    --i-type fasta 2>&1
grep -q 'HiSeq2500_150bp' "${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log
AM_EXEC_SHOULD_FAIL \
    --i-file "${OUT_DIR}"/devnull \
    --lc "pe" \
    --i-parser memory \
    --builtin_qual_file MiniSeq_TruSeq_50bp \
    --i-fcov 0 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 200 \
    --i-type fasta 2>&1
grep -q 'MiniSeq_TruSeq_50bp is not a valid paired-end profile' "${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log
AM_EXEC_SHOULD_FAIL \
    --i-file "${OUT_DIR}"/devnull \
    --lc "pe" \
    --i-parser memory \
    --qual_file_1 "${OUT_DIR}"/devnull --qual_file_2 "${OUT_DIR}"/devnull \
    --i-fcov 0 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 200 \
    --i-type fasta 2>&1
grep -q 'Profile empty' "${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log
AM_EXEC_SHOULD_FAIL \
    --i-file "${OUT_DIR}"/devnull \
    --lc "pe" \
    --i-parser memory \
    --qual_file_1 data/Illumina_profiles/Emp36R1.txt \
    --i-fcov 0 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 200 \
    --i-type fasta 2>&1
grep -q 'An input file path for --qual_file_2 must be specified.' "${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log
AM_EXEC \
    --i-file "${OUT_DIR}"/devnull \
    --lc "pe" \
    --i-parser memory \
    --qual_file_1 data/Illumina_profiles/Emp36R1.txt \
    --qual_file_2 data/Illumina_profiles/Emp36R2.txt \
    --i-fcov 0 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 200 \
    --i-type fasta
grep -q 'Read quality profile size for R1: 36' "${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log

unset KEEPLOG
rm -fr "${OUT_DIR}"/log_*.d "${OUT_DIR}"/am_exec_*.log "${OUT_DIR}"/devnull
assert_cleandir
