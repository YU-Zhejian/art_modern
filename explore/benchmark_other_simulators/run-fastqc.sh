#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
# shellcheck disable=SC2155
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

export OUT_DIR="$(pwd)"/data_out
export ART_MODERN_THREADS=20
export FCOV=3

mkdir -p "${OUT_DIR}"

bin/wgsim \
    -1 100 -2 100 -N 185412 -d 300 -s 20 -r 0 \
    data/yeast.fa \
    "${OUT_DIR}"/yeast_wgsim_1.fq "${OUT_DIR}"/yeast_wgsim_2.fq &

bin/dwgsim \
    -1 100 -2 100 -C "${FCOV}" -d 300 -s 20 -o 1 -r 0 -y 0 \
    data/yeast.fa \
    "${OUT_DIR}"/yeast_dwgsim &

bin/art_original \
    --in data/yeast.fa --out "${OUT_DIR}"/yeast_art_ \
    --qprof1 data/e_coli_art_R1.txt \
    --qprof2 data/e_coli_art_R2.txt \
    --fcov "${FCOV}" --len 100 --mflen 300 --sdev 20 --noALN --paired &

bin/pirs simulate -A dist -m 300 -l 100 -x "${FCOV}" -v 20 -t 20 \
    -B data/e_coli_HiSeq2K_pirs_bcm.count.matrix \
    -I data/e_coli_pirs_indelstat.InDel.matrix \
    --no-gc-bias \
    -o "${OUT_DIR}"/pirs \
    -c text \
    data/yeast.fa
mv "${OUT_DIR}"/pirs/Sim_100_300_1.fq "${OUT_DIR}"/pirs_1.fq
mv "${OUT_DIR}"/pirs/Sim_100_300_2.fq "${OUT_DIR}"/pirs_2.fq
wait
opt/art_modern_build/art_modern \
    --mode wgs --lc pe \
    --i-file data/yeast.fa --i-fcov "${FCOV}" --read_len 100 \
    --o-fastq "${OUT_DIR}"/yeast_art_modern_wgs_memory.fastq \
    --qual_file_1 data/e_coli_art_R1.txt \
    --qual_file_2 data/e_coli_art_R2.txt \
    --pe_frag_dist_mean 300 --pe_frag_dist_std_dev 20 --parallel "${ART_MODERN_THREADS}"
seqtk seq -1 "${OUT_DIR}"/yeast_art_modern_wgs_memory.fastq >"${OUT_DIR}"/yeast_art_modern_wgs_memory_1.fastq
seqtk seq -2 "${OUT_DIR}"/yeast_art_modern_wgs_memory.fastq >"${OUT_DIR}"/yeast_art_modern_wgs_memory_2.fastq
rm -f "${OUT_DIR}"/yeast_art_modern_wgs_memory.fastq

fastqc \
    --extract \
    --svg \
    --threads 20 \
    "${OUT_DIR}"/*.fastq \
    "${OUT_DIR}"/*.fq &
# rm -fr "${OUT_DIR}"/*.fastq "${OUT_DIR}"/*.fq

mkdir -p profiles
python get_profile.py data/e_coli_SCNR0028307_p33_1.fastq.gz profiles/real_1.profile &
python get_profile.py data/e_coli_SCNR0028307_p33_2.fastq.gz profiles/real_2.profile &
python get_profile.py data_out/yeast_art_modern_wgs_memory_1.fastq profiles/art_modern_1.profile &
python get_profile.py data_out/yeast_art_modern_wgs_memory_2.fastq profiles/art_modern_2.profile &
python get_profile.py data_out/yeast_art_1.fq profiles/art_original_1.profile &
python get_profile.py data_out/yeast_art_2.fq profiles/art_original_2.profile &
python get_profile.py data_out/pirs_1.fq profiles/pirs_1.profile &
python get_profile.py data_out/pirs_2.fq profiles/pirs_2.profile &
python get_profile.py data_out/yeast_dwgsim.bwa.read1.fastq profiles/dwgsim_1.profile &
python get_profile.py data_out/yeast_dwgsim.bwa.read2.fastq profiles/dwgsim_2.profile &
python get_profile.py data_out/yeast_wgsim_1.fq profiles/wgsim_1.profile &
python get_profile.py data_out/yeast_wgsim_2.fq profiles/wgsim_2.profile &
wait
python cmp_profile.py >correlation.tsv
