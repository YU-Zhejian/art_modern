# shellcheck shell=bash

parser=memory
coverage=pbsim3
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file data/raw_data/ce11.mRNA_head.pbsim3.transcript \
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
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}" "data/raw_data/ce11.mRNA_head.fa"
                python sh.d/test-small.sh.d/test_sam.py \
                    "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq \
                    data/raw_data/ce11.mRNA_head.pbsim3.transcript \
                    data/raw_data/ce11.mRNA_head.pbsim3.transcript \
                    PBSIM3_TRANSCRIPT
done

parser=stream
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file data/raw_data/ce11.mRNA_head.pbsim3.transcript \
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
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl "data/raw_data/ce11.mRNA_head.fa"
                    python sh.d/test-small.sh.d/test_sam.py \
                        "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".fq \
                        data/raw_data/ce11.mRNA_head.pbsim3.transcript \
                        data/raw_data/ce11.mRNA_head.pbsim3.transcript \
                        PBSIM3_TRANSCRIPT
done
rm -fr "${OUT_DIR}"/test_small_??_template_"${parser}"_"${coverage}".hl.sam
