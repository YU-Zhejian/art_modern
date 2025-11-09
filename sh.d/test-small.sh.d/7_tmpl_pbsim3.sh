# shellcheck shell=bash

parser=memory
coverage=pbsim3
for lc in se pe mp; do
    AM_EXEC \
        --i-file "${MRNA_PBSIM3_TRANSCRIPT}" \
        --i-type pbsim3_transcripts \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
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
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        PBSIM3_TRANSCRIPT \
        IS_TEMPLATE
    rm -fr "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq
    assert_cleandir
done

parser=stream
for lc in se pe mp; do
    AM_EXEC \
        --i-file "${MRNA_PBSIM3_TRANSCRIPT}" \
        --i-type pbsim3_transcripts \
        --i-batch_size 100 \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
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
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        PBSIM3_TRANSCRIPT \
        IS_TEMPLATE
    rm -fr \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq \
        "${OUT_DIR}"test_small_"${lc}"_template_"${parser}"_"${coverage}".hl.sam
    assert_cleandir
done
rm -fr "${OUT_DIR}"/test_small_??_template_"${parser}"_"${coverage}".hl.sam

# We would now test whether our separated --read_len_1 and --read_len_2 is working.
