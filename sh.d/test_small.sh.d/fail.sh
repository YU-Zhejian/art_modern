# shell=bash
FCOV=10
parser=memory
for lc in se pe mp; do
    for mode in wgs trans template; do
        "${ART}" \
            --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt \
            --qual_file_2 data/Illumina_profiles/HiSeq2500L125R2.txt \
            --i-file "${ERR_FA}" \
            --read_len 125 \
            --mode "${mode}" \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov "${FCOV}" \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-fastq "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq \
            --pe_frag_dist_std_dev 20 \
            --pe_frag_dist_mean 200
        cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
        rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq
    done
done

unset FCOV
