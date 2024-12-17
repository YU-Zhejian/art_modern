#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

export OUT_DIR=/tmp/data_out
export ART_MODERN_THREADS=20

mkdir -p "${OUT_DIR}"

function run() {
    /bin/time -a -o time.tsv -f "${1}"'\t%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c' "${@:2}"
}
function recreate_data_out() {
    rm -rf "${OUT_DIR:?}"/*
    mkdir -p "${OUT_DIR}"
}

printf 'TEST_CASE\tWALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\n' >time.tsv
for i in {1..3}; do
    echo "Run ${i}"

    run wgsim-genome bin/wgsim \
        -1 150 -2 150 -N 3409740 -d 300 -s 20 -r 0 \
        data/ce11.fa \
        "${OUT_DIR}"/ce11_wgsim_1.fq "${OUT_DIR}"/ce11_wgsim_2.fq
    recreate_data_out

    run dwgsim-genome bin/dwgsim \
        -1 150 -2 150 -C 10 -d 300 -s 20 -o 2 -r 0 -y 0 \
        data/ce11.fa \
        "${OUT_DIR}"/ce11_dwgsim
    recreate_data_out

    run art_original-genome bin/art_original \
        --in data/ce11.fa --out "${OUT_DIR}"/ce11_art_ \
        -f 10 --len 150 --mflen 300 --sdev 20 --noALN --paired --seqSys HS25
    recreate_data_out

    run art_modern-genome opt/art_modern_build/art_modern \
        --mode wgs --lc pe \
        --i-file data/ce11.fa --i-fcov 10 --read_len 150 \
        --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
        --qual_file_1 ../../data/Illumina_profiles/HiSeq2500L150R1filter.txt \
        --qual_file_2 ../../data/Illumina_profiles/HiSeq2500L150R2filter.txt \
        --pe_frag_dist_mean 300 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"
    recreate_data_out

    run wgsim-transcriptome bin/wgsim \
        -1 150 -2 150 -N 13393233 -d 300 -s 20 -r 0 \
        data/hg38_long_mrna.fa \
        "${OUT_DIR}"/hg38_long_mrna_wgsim_1.fq "${OUT_DIR}"/hg38_long_mrna_wgsim_2.fq
    recreate_data_out

    run art_original-transcriptome bin/art_original \
        --in data/hg38_long_mrna.fa --out "${OUT_DIR}"/hg38_long_mrna_art_ \
        -f 4 --len 150 --mflen 300 --sdev 20 --noALN --paired --seqSys HS25
    recreate_data_out

    run art_modern-transcriptome opt/art_modern_build/art_modern \
        --mode trans --lc pe \
        --i-file data/hg38_long_mrna.fa --i-fcov 4 --read_len 150 --i-parser stream \
        --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
        --qual_file_1 ../../data/Illumina_profiles/HiSeq2500L150R1filter.txt \
        --qual_file_2 ../../data/Illumina_profiles/HiSeq2500L150R2filter.txt \
        --pe_frag_dist_mean 300 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
        --i-batch_size 1024
    recreate_data_out
done
rm -fr "${OUT_DIR}"
