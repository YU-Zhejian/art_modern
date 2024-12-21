#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
# shellcheck disable=SC2155
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

export OUT_DIR="$(pwd)"/data_out
export FCOV=3

rm -fr "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

bin/wgsim \
    -1 100 -2 100 -N "$(echo "${FCOV} * $(wc -c <data/yeast.fa) / 200" | bc -q)" -d 500 -s 20 -r 0 \
    data/yeast.fa \
    "${OUT_DIR}"/pe100_wgsim_1.fq "${OUT_DIR}"/pe100_wgsim_2.fq &
bin/wgsim \
    -1 300 -2 300 -N "$(echo "${FCOV} * $(wc -c <data/yeast.fa) / 600" | bc -q)" -d 500 -s 20 -r 0 \
    data/yeast.fa \
    "${OUT_DIR}"/pe300_wgsim_1.fq "${OUT_DIR}"/pe300_wgsim_2.fq &

bin/dwgsim \
    -1 100 -2 100 -C "${FCOV}" -d 500 -s 20 -o 1 -r 0 -y 0 \
    data/yeast.fa \
    "${OUT_DIR}"/pe100_dwgsim &
bin/dwgsim \
    -1 300 -2 300 -C "${FCOV}" -d 500 -s 20 -o 1 -r 0 -y 0 \
    data/yeast.fa \
    "${OUT_DIR}"/pe300_dwgsim &

bin/art_original \
    --in data/yeast.fa --out "${OUT_DIR}"/pe100_art_ \
    --qprof1 data/e_coli_HiSeq2K_art_R1.txt \
    --qprof2 data/e_coli_HiSeq2K_art_R2.txt \
    --fcov "${FCOV}" --len 100 --mflen 500 --sdev 20 --noALN --paired &
bin/art_original \
    --in data/yeast.fa --out "${OUT_DIR}"/pe300_art_ \
    --qprof1 data/soybean_HiSeq2500_art_R1.txt \
    --qprof2 data/soybean_HiSeq2500_art_R2.txt \
    --fcov "${FCOV}" --len 300 --mflen 500 --sdev 20 --noALN --paired &

bin/pirs simulate -A dist -m 500 -l 100 -x "${FCOV}" -v 20 -t 20 \
    -B data/e_coli_HiSeq2K_pirs_bcm.count.matrix \
    -I data/e_coli_HiSeq2K_pirs_indelstat.InDel.matrix \
    --no-gc-bias \
    -o "${OUT_DIR}"/pirs \
    -c text \
    data/yeast.fa &
bin/pirs simulate -A dist -m 500 -l 300 -x "${FCOV}" -v 20 -t 20 \
    -B data/e_coli_HiSeq2K_pirs_bcm.count.matrix \
    -I data/e_coli_HiSeq2K_pirs_indelstat.InDel.matrix \
    --no-gc-bias \
    -o "${OUT_DIR}"/pirs \
    -c text \
    data/yeast.fa &

opt/art_modern_build/art_modern --lc pe \
    --i-file data/yeast.fa --i-fcov "${FCOV}" --read_len 100 \
    --o-fastq "${OUT_DIR}"/pe100_art_modern.fq \
    --qual_file_1 data/e_coli_HiSeq2K_art_R1.txt \
    --qual_file_2 data/e_coli_HiSeq2K_art_R2.txt \
    --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 &
opt/art_modern_build/art_modern --lc pe \
    --i-file data/yeast.fa --i-fcov "${FCOV}" --read_len 300 \
    --o-fastq "${OUT_DIR}"/pe300_art_modern.fq \
    --qual_file_1 data/soybean_HiSeq2500_art_R1.txt \
    --qual_file_2 data/soybean_HiSeq2500_art_R2.txt \
    --pe_frag_dist_mean 500 --pe_frag_dist_std_dev 20 &

wait
for i in 1 2; do
    for rlen in 100 300; do
        mv "${OUT_DIR}"/pirs/Sim_"${rlen}"_500_"${i}".fq "${OUT_DIR}"/pe"${rlen}"_pirs_"${i}".fq

        mv "${OUT_DIR}"/pe"${rlen}"_dwgsim.bwa.read"${i}".fastq "${OUT_DIR}"/pe"${rlen}"_dwgsim_"${i}".fq

        seqtk seq -"${i}" "${OUT_DIR}"/pe"${rlen}"_art_modern.fq >"${OUT_DIR}"/pe"${rlen}"_art_modern_"${i}".fq
    done
done

rm -fr \
    "${OUT_DIR}"/pe100_art_modern.fq \
    "${OUT_DIR}"/pe300_art_modern.fq \
    "${OUT_DIR}"/pirs \
    "${OUT_DIR}"/*.mutations.*

mkdir -p "${OUT_DIR}"/fastqc
fastqc \
    --outdir "${OUT_DIR}"/fastqc \
    --extract \
    --svg \
    --threads 20 \
    "${OUT_DIR}"/*.fq &

mkdir -p profiles
for rlen in 100 300; do
    for software in art_modern art pirs wgsim dwgsim; do
        for i in 1 2; do
            python get_profile.py \
                "${OUT_DIR}"/pe"${rlen}"_"${software}"_"${i}".fq \
                profiles/pe"${rlen}"_"${software}"_"${i}".profile "${rlen}" &
        done
    done
done
wait

python cmp_profile.py >correlation.tsv
Rscript plot_corr.R
