# shellcheck shell=bash

AM_EXEC \
    --builtin_qual_file HiSeq2500_150bp \
    --sep_flag \
    --i-file data/raw_data/ce11_chr1.fa \
    --read_len 150 \
    --mode wgs \
    --lc se \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq "${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq \
    --o-fasta "${OUT_DIR}"/test_small_se_wgs_memory_sep.fasta \
    --o-pwa "${OUT_DIR}"/test_small_se_wgs_memory_sep.pwa \
    --o-sam "${OUT_DIR}"/test_small_se_wgs_memory_sep.sam \
    --o-hl_sam "${OUT_DIR}"/test_small_se_wgs_memory_sep.hl.bam \
    --o-hl_sam-num_threads 2 \
    --o-hl_sam-compress_level u \
    --o-hl_sam-write_bam
merge_file "${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq
merge_file "${OUT_DIR}"/test_small_se_wgs_memory_sep.sam
merge_file "${OUT_DIR}"/test_small_se_wgs_memory_sep.fasta

if type -p fastqc &>/dev/null && type -p x-www-browser &>/dev/null && [ ! "${NO_FASTQC:-}" = "1" ]; then
    fastqc -t "${PARALLEL}" "${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq
    # Open the browser and ignore what's happening afterwards
    {
        x-www-browser "${OUT_DIR}"/test_small_se_wgs_memory_sep_fastqc.html &>/dev/null || true
    } &
    sleep 3
fi

samtools fastq "${OUT_DIR}"/test_small_se_wgs_memory_sep.sam | seqkit sort >"${OUT_DIR}"/test_small_se_wgs_memory_sep.sam.fq
seqkit sort <"${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq >"${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq.srt.fq
cmp \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq.srt.fq \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep.sam.fq
seqkit sort --line-width 0 <"${OUT_DIR}"/test_small_se_wgs_memory_sep.fasta >"${OUT_DIR}"/test_small_se_wgs_memory_sep.fasta.srt.fa
seqtk seq -A <"${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq.srt.fq >"${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq.srt.fa
cmp \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq.srt.fa \
    "${OUT_DIR}"/test_small_se_wgs_memory_sep.fasta.srt.fa

rm -fr "${OUT_DIR}"/test_small_se_wgs_memory_sep.* "${OUT_DIR}"/test_small_se_wgs_memory_sep_fastqc.*
assert_cleandir

AM_EXEC \
    --builtin_qual_file HiSeq2500_150bp \
    --sep_flag \
    --i-file "${LAMBDA_PHAGE}" \
    --read_len 150 \
    --mode wgs \
    --lc pe \
    --i-parser memory \
    --i-fcov 5 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq "${OUT_DIR}"/test_small_pe_wgs_memory_sep.fastq \
    --o-fasta "${OUT_DIR}"/test_small_pe_wgs_memory_sep.fasta \
    --o-pwa "${OUT_DIR}"/test_small_pe_wgs_memory_sep.pwa \
    --o-sam "${OUT_DIR}"/test_small_pe_wgs_memory_sep.sam \
    --o-hl_sam "${OUT_DIR}"/test_small_pe_wgs_memory_sep.hl.bam \
    --o-hl_sam-num_threads 2 \
    --o-hl_sam-compress_level u \
    --o-hl_sam-write_bam \
    --pe_frag_dist_std_dev 20 \
    --pe_frag_dist_mean 500
merge_file "${OUT_DIR}"/test_small_pe_wgs_memory_sep.fastq
merge_file "${OUT_DIR}"/test_small_pe_wgs_memory_sep.sam
merge_file "${OUT_DIR}"/test_small_pe_wgs_memory_sep.fasta

samtools fastq \
    -1 "${OUT_DIR}"/test_small_pe_wgs_memory_sep.sam.1.fq \
    -2 "${OUT_DIR}"/test_small_pe_wgs_memory_sep.sam.2.fq \
    -N \
    "${OUT_DIR}"/test_small_pe_wgs_memory_sep.sam
for fn in \
    "${OUT_DIR}"/test_small_pe_wgs_memory_sep.sam.1.fq \
    "${OUT_DIR}"/test_small_pe_wgs_memory_sep.sam.2.fq; do
    seqkit sort <"$fn" >"$fn".srt.fq
    mv "$fn".srt.fq "$fn"
done
seqtk seq -1 "${OUT_DIR}"/test_small_pe_wgs_memory_sep.fastq |
    seqkit sort >"${OUT_DIR}"/test_small_pe_wgs_memory_sep.fastq.1.fq
seqtk seq -2 "${OUT_DIR}"/test_small_pe_wgs_memory_sep.fastq |
    seqkit sort >"${OUT_DIR}"/test_small_pe_wgs_memory_sep.fastq.2.fq

for i in 1 2; do
    cmp \
        "${OUT_DIR}"/test_small_pe_wgs_memory_sep.fastq."${i}".fq \
        "${OUT_DIR}"/test_small_pe_wgs_memory_sep.sam."${i}".fq
done

rm -fr "${OUT_DIR}"/test_small_pe_wgs_memory_sep.*
assert_cleandir
# I suppose the FASTA format does not need to be tested here.

if [ "${FORMAT_ONLY:-}" = "1" ]; then
    exit 0
fi
