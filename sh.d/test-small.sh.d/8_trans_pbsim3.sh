# shellcheck shell=bash

# Trans mode with stranded/strandless coverage
coverage=pbsim3
parser=memory
for lc in se pe mp; do
    AM_EXEC \
        --i-file "${MRNA_PBSIM3_TRANSCRIPT}" \
        --i-type pbsim3_transcripts \
        --mode trans \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".sam \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".fq \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 500
    merge_file "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".sam
    merge_file "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".fq
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}" "${MRNA_HEAD}"
    python sh.d/test-small.sh.d/validate_cov.py \
        "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".fq \
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        PBSIM3_TRANSCRIPT \
        NOT_TEMPLATE
    rm -fr "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".fq
    assert_cleandir
done
