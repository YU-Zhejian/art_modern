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
    echo "${1}"
    /bin/time -a -o time.tsv -f "${1}"'\t%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c' "${@:2}" &>"${OUT_DIR}"/"${1}".log || return 1
    rm -rf "${OUT_DIR:?}"/*
    mkdir -p "${OUT_DIR}"
}

printf 'TEST_CASE\tWALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\n' >time.tsv
for i in {1..3}; do
    echo "Run ${i}"

    run pirs-genome-pe100 bin/pirs simulate -A dist -m 500 -l 100 -x 10 -v 20 -t 20 \
        -B data/e_coli_HiSeq2K_pirs_bcm.count.matrix \
        -I data/e_coli_HiSeq2K_pirs_indelstat.InDel.matrix \
        --no-gc-bias \
        -o "${OUT_DIR}"/Illumina \
        -c text \
        data/ce11.fa

    run dwgsim-genome-pe100 bin/dwgsim \
        -1 100 -2 100 -C 10 -d 500 -s 20 -o 2 -r 0 -y 0 \
        data/ce11.fa \
        "${OUT_DIR}"/ce11_dwgsim

    run wgsim-genome-pe100 bin/wgsim \
        -1 100 -2 100 -N "$(echo "10 * $(wc -c <data/ce11.fa) / 200" | bc -q)" -d 500 -s 20 -r 0 \
        data/ce11.fa \
        "${OUT_DIR}"/ce11_wgsim_1.fq "${OUT_DIR}"/ce11_wgsim_2.fq

    run art_original-genome-pe100 bin/art_original \
        --in data/ce11.fa --out "${OUT_DIR}"/ce11_art_ \
        --qprof1 data/e_coli_HiSeq2K_art_R1.txt \
        --qprof2 data/e_coli_HiSeq2K_art_R2.txt \
        -f 10 --len 100 --mflen 500 --sdev 20 --noALN --paired

    run art_modern-genome-pe100 opt/art_modern_build/art_modern \
        --mode wgs --lc pe \
        --i-file data/ce11.fa --i-fcov 10 --read_len 100 \
        --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
        --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
        --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"

    run art_modern_gcc-genome-pe100 opt/art_modern_gcc_build/art_modern \
        --mode wgs --lc pe \
        --i-file data/ce11.fa --i-fcov 10 --read_len 100 \
        --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
        --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
        --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"

    run art_modern_prev-genome-pe100 opt/art_modern_prev_ver_build/art_modern \
        --mode wgs --lc pe \
        --i-file data/ce11.fa --i-fcov 10 --read_len 100 \
        --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
        --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
        --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"

    #    run pirs-transcriptome-pe100 bin/pirs simulate -A dist -m 500 -l 100 -x 4 -v 20 -t 20 \
    #        -B data/e_coli_HiSeq2K_pirs_bcm.count.matrix \
    #        -I data/e_coli_HiSeq2K_pirs_indelstat.InDel.matrix \
    #        --no-gc-bias \
    #        -o "${OUT_DIR}"/Illumina \
    #        -c text \
    #        data/hg38_long_mrna.fa

    run dwgsim-transcriptome-pe100 bin/dwgsim \
        -1 100 -2 100 -C 4 -d 500 -s 20 -o 2 -r 0 -y 0 \
        data/hg38_long_mrna.fa \
        "${OUT_DIR}"/hg38_long_mrna_dwgsim

    run wgsim-transcriptome-pe100 bin/wgsim \
        -1 100 -2 100 -N "$(echo "4 * $(wc -c <data/hg38_long_mrna.fa) / 200" | bc -q)" -d 500 -s 20 -r 0 \
        data/hg38_long_mrna.fa \
        "${OUT_DIR}"/hg38_long_mrna_wgsim_1.fq "${OUT_DIR}"/hg38_long_mrna_wgsim_2.fq

    run art_original-transcriptome-pe100 bin/art_original \
        --in data/hg38_long_mrna.fa --out "${OUT_DIR}"/hg38_long_mrna_art_ \
        --qprof1 data/e_coli_HiSeq2K_art_R1.txt \
        --qprof2 data/e_coli_HiSeq2K_art_R2.txt \
        -f 4 --len 100 --mflen 500 --sdev 20 --noALN --paired

    run art_modern-transcriptome-pe100 opt/art_modern_build/art_modern \
        --mode trans --lc pe \
        --i-file data/hg38_long_mrna.fa --i-fcov 4 --read_len 100 --i-parser stream \
        --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
        --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
        --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
        --i-batch_size 1024

    run art_modern_gcc-transcriptome-pe100 opt/art_modern_gcc_build/art_modern \
        --mode trans --lc pe \
        --i-file data/hg38_long_mrna.fa --i-fcov 4 --read_len 100 --i-parser stream \
        --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
        --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
        --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
        --i-batch_size 1024

    run art_modern_prev-transcriptome-pe100 opt/art_modern_prev_ver_build/art_modern \
        --mode trans --lc pe \
        --i-file data/hg38_long_mrna.fa --i-fcov 4 --read_len 100 --i-parser stream \
        --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
        --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
        --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
        --i-batch_size 1024

    # Now run the same tests with PE300 coverage

    run pirs-genome-pe300 bin/pirs simulate -A dist -m 500 -l 300 -x 10 -v 20 -t 20 \
        -B data/soybean_HiSeq2500_pirs_bcm.count.matrix \
        -I data/soybean_HiSeq2500_pirs_indelstat.InDel.matrix \
        --no-gc-bias \
        -o "${OUT_DIR}"/Illumina \
        -c text \
        data/ce11.fa

    run dwgsim-genome-pe300 bin/dwgsim \
        -1 300 -2 300 -C 10 -d 500 -s 20 -o 2 -r 0 -y 0 \
        data/ce11.fa \
        "${OUT_DIR}"/ce11_dwgsim

    run wgsim-genome-pe300 bin/wgsim \
        -1 300 -2 300 -N "$(echo "10 * $(wc -c <data/ce11.fa) / 200" | bc -q)" -d 500 -s 20 -r 0 \
        data/ce11.fa \
        "${OUT_DIR}"/ce11_wgsim_1.fq "${OUT_DIR}"/ce11_wgsim_2.fq

    run art_original-genome-pe300 bin/art_original \
        --in data/ce11.fa --out "${OUT_DIR}"/ce11_art_ \
        --qprof1 data/soybean_HiSeq2500_art_R1.txt \
        --qprof2 data/soybean_HiSeq2500_art_R2.txt \
        -f 10 --len 300 --mflen 500 --sdev 20 --noALN --paired

    run art_modern-genome-pe300 opt/art_modern_build/art_modern \
        --mode wgs --lc pe \
        --i-file data/ce11.fa --i-fcov 10 --read_len 300 \
        --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
        --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
        --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"

    run art_modern_gcc-genome-pe300 opt/art_modern_gcc_build/art_modern \
        --mode wgs --lc pe \
        --i-file data/ce11.fa --i-fcov 10 --read_len 300 \
        --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
        --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
        --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"

    run art_modern_prev-genome-pe300 opt/art_modern_prev_ver_build/art_modern \
        --mode wgs --lc pe \
        --i-file data/ce11.fa --i-fcov 10 --read_len 300 \
        --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
        --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
        --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"

    #    run pirs-transcriptome-pe300 bin/pirs simulate -A dist -m 500 -l 300 -x 4 -v 20 -t 20 \
    #        -B data/soybean_HiSeq2500_pirs_bcm.count.matrix \
    #        -I data/soybean_HiSeq2500_pirs_indelstat.InDel.matrix \
    #        --no-gc-bias \
    #        -o "${OUT_DIR}"/Illumina \
    #        -c text \
    #        data/hg38_long_mrna.fa

    run dwgsim-transcriptome-pe300 bin/dwgsim \
        -1 300 -2 300 -C 4 -d 500 -s 20 -o 2 -r 0 -y 0 \
        data/hg38_long_mrna.fa \
        "${OUT_DIR}"/hg38_long_mrna_dwgsim

    run wgsim-transcriptome-pe300 bin/wgsim \
        -1 300 -2 300 -N "$(echo "4 * $(wc -c <data/hg38_long_mrna.fa) / 200" | bc -q)" -d 500 -s 20 -r 0 \
        data/hg38_long_mrna.fa \
        "${OUT_DIR}"/hg38_long_mrna_wgsim_1.fq "${OUT_DIR}"/hg38_long_mrna_wgsim_2.fq

    run art_original-transcriptome-pe300 bin/art_original \
        --in data/hg38_long_mrna.fa --out "${OUT_DIR}"/hg38_long_mrna_art_ \
        --qprof1 data/soybean_HiSeq2500_art_R1.txt \
        --qprof2 data/soybean_HiSeq2500_art_R2.txt \
        -f 4 --len 300 --mflen 500 --sdev 20 --noALN --paired

    run art_modern-transcriptome-pe300 opt/art_modern_build/art_modern \
        --mode trans --lc pe \
        --i-file data/hg38_long_mrna.fa --i-fcov 4 --read_len 300 --i-parser stream \
        --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
        --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
        --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
        --i-batch_size 1024

    run art_modern_gcc-transcriptome-pe300 opt/art_modern_gcc_build/art_modern \
        --mode trans --lc pe \
        --i-file data/hg38_long_mrna.fa --i-fcov 4 --read_len 300 --i-parser stream \
        --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
        --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
        --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
        --i-batch_size 1024

    run art_modern_prev-transcriptome-pe300 opt/art_modern_prev_ver_build/art_modern \
        --mode trans --lc pe \
        --i-file data/hg38_long_mrna.fa --i-fcov 4 --read_len 300 --i-parser stream \
        --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
        --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
        --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
        --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}" \
        --i-batch_size 1024
done
rm -fr "${OUT_DIR}"
