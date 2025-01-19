# shell=bash
# Trans mode with constant coverage
parser=memory
FCOV=10
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file "${MRNA_HEAD}" \
        --read_len 125 \
        --mode trans \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}".sam \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 500
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_trans_"${parser}" "${MRNA_HEAD}"
done
# No need to test stream FASTA parser
unset FCOV
