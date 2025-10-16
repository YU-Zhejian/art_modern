# shellcheck shell=bash

# Trans mode with stranded/strandless coverage
parser=memory
for coverage in stranded strandless; do
    for lc in se pe mp; do
        "${ART}" \
            --builtin_qual_file HiSeq2500_125bp \
            --i-file "${MRNA_HEAD}" \
            --read_len 125 \
            --mode trans \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov data/raw_data/ce11.mRNA_head.cov_"${coverage}".tsv \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-sam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".sam \
            --o-fastq "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".fq \
            --pe_frag_dist_std_dev 20 \
            --pe_frag_dist_mean 500
        sam2bam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}" "${MRNA_HEAD}"
        python sh.d/test-small.sh.d/validate_cov.py \
            "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}"_"${coverage}".fq \
            "${MRNA_HEAD}" \
            data/raw_data/ce11.mRNA_head.cov_"${coverage}".tsv \
            COV_TSV \
            NOT_TEMPLATE
    done
done
# No need to test stream FASTA parser
