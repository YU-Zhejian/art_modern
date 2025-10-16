# shellcheck shell=bash

parser=memory
coverage=pbsim3
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file "${MRNA_PBSIM3_TRANSCRIPT}" \
        --read_len 125 \
        --i-type pbsim3_transcripts \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".sam \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}" "${MRNA_HEAD}"
    python sh.d/test-small.sh.d/validate_cov.py \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq \
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        PBSIM3_TRANSCRIPT \
        IS_TEMPLATE
done

parser=stream
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file "${MRNA_PBSIM3_TRANSCRIPT}" \
        --read_len 125 \
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
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl "${MRNA_HEAD}"
    python sh.d/test-small.sh.d/validate_cov.py \
        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq \
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        "${MRNA_PBSIM3_TRANSCRIPT}" \
        PBSIM3_TRANSCRIPT \
        IS_TEMPLATE
done
rm -fr "${OUT_DIR}"/test_small_??_template_"${parser}"_"${coverage}".hl.sam
