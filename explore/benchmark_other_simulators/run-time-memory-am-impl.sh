#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

export OUT_DIR=data_out_am_impl
export ART_MODERN_THREADS=20

mkdir -p "${OUT_DIR}"

function run() {
    echo "${1}"
    /bin/time -a -o time.tsv -f "${1}"'\t%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c' "${@:2}" &>"${OUT_DIR}"/"${1}".log || return 1
    # rm -rf "${OUT_DIR:?}"/*
    mkdir -p "${OUT_DIR}"
}

printf 'TEST_CASE\tWALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\n' >time.tsv
for i in {1..3}; do
    echo "Run ${i}"
    for name in art_modern art_modern_stlmap art_modern_gcc art_modern_prev_ver; do

        run "${name}"-genome-pe100 opt/"${name}"_build/art_modern \
            --mode wgs --lc pe \
            --i-file data/ce11.fa --i-fcov 100 --read_len 100 \
            --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
            --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
            --o-fastq /dev/null \
            --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"

        run "${name}"-transcriptome-pe100 opt/"${name}"_build/art_modern \
            --mode trans --lc pe \
            --i-file data/hg38_long_mrna.fa --i-fcov 40 --read_len 100 --i-parser stream \
            --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
            --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
            --o-fastq /dev/null \
            --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
            --i-batch_size 1024

        # Now run the same tests with PE300 coverage

        run "${name}"-genome-pe300 opt/"${name}"_build/art_modern \
            --mode wgs --lc pe \
            --i-file data/ce11.fa --i-fcov 100 --read_len 300 \
            --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
            --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
            --o-fastq /dev/null \
            --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"

        run "${name}"-transcriptome-pe300 opt/"${name}"_build/art_modern \
            --mode trans --lc pe \
            --i-file data/hg38_long_mrna.fa --i-fcov 40 --read_len 300 --i-parser stream \
            --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
            --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
            --o-fastq /dev/null \
            --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
            --i-batch_size 1024
    done
done
rm -fr "${OUT_DIR}"
