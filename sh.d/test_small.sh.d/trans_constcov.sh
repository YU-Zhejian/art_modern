# shell=bash
# Trans mode with constant coverage
parser=memory
for lc in se pe mp; do
    build/art_modern \
        --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
        --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
        --i-file raw_data/ce11.mRNA_head.fa \
        --read_len 125 \
        --mode trans \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam tmp/test_small_"${lc}"_trans_"${parser}".sam \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 500
    sam2bam tmp/test_small_"${lc}"_trans_"${parser}" raw_data/ce11.mRNA_head.fa
done
# No need to test stream FASTA parser
