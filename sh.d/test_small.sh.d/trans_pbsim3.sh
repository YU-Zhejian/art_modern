# shell=bash
# Trans mode with stranded/strandless coverage
coverage=pbsim3
parser=memory
for lc in se pe mp; do
    build/art_modern \
        --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
        --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
        --i-file raw_data/ce11.mRNA_head.pbsim3.transcript \
        --read_len 125 \
        --i-type pbsim3_template \
        --mode trans \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam tmp/test_small_"${lc}"_trans_"${parser}"_"${coverage}".sam \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 500
    sam2bam tmp/test_small_"${lc}"_trans_"${parser}"_"${coverage}" raw_data/ce11.mRNA_head.fa
done
