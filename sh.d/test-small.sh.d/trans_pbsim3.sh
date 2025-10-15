# shellcheck shell=bash

# Trans mode with stranded/strandless coverage
coverage=pbsim3
parser=memory
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file data/raw_data/ce11.mRNA_head.pbsim3.transcript \
        --read_len 125 \
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
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}" "${MRNA_HEAD}"
            python sh.d/test-small.sh.d/test_sam.py \
                "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".fq \
                data/raw_data/ce11.mRNA_head.pbsim3.transcript \
                data/raw_data/ce11.mRNA_head.pbsim3.transcript \
                PBSIM3_TRANSCRIPT
done
