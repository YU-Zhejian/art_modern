# shellcheck shell=bash
FCOV=10
parser=memory
err_file="data/raw_data/err.fa"
for lc in se pe mp; do
    for mode in wgs trans template; do
        "${ART}" \
            --builtin_qual_file HiSeq2500_125bp \
            --i-file "${err_file}" \
            --read_len 125 \
            --mode "${mode}" \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov "${FCOV}" \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-fastq "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq \
            --o-fasta "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fasta \
            --o-pwa "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.pwa \
            --o-sam "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.sam \
            --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.hl.bam \
            --o-hl_sam-num_threads 2 \
            --o-hl_sam-compress_level u \
            --o-hl_sam-write_bam \
            --pe_frag_dist_std_dev 20 \
            --pe_frag_dist_mean 200 \
            --i-type fasta
        cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
        rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
    done
done

parser=htslib
mode=wgs
for lc in se pe mp; do
    "${ART}" \
        --builtin_qual_file HiSeq2500_125bp \
        --i-file "${err_file}" \
        --read_len 125 \
        --mode "${mode}" \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --i-fcov "${FCOV}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-fastq "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq \
        --o-fasta "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fasta \
        --o-pwa "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.pwa \
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.sam \
        --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.hl.bam \
        --o-hl_sam-num_threads 2 \
        --o-hl_sam-compress_level u \
        --o-hl_sam-write_bam \
        --pe_frag_dist_std_dev 20 \
        --pe_frag_dist_mean 200 \
        --i-type fasta
    cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
    rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
done

parser=stream
for lc in se pe mp; do
    for mode in trans template; do
        "${ART}" \
            --builtin_qual_file HiSeq2500_125bp \
            --i-file "${err_file}" \
            --read_len 125 \
            --mode "${mode}" \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov "${FCOV}" \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-fastq "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq \
            --o-fasta "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fasta \
            --o-pwa "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.pwa \
            --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.hl.bam \
            --o-hl_sam-num_threads 2 \
            --o-hl_sam-compress_level u \
            --o-hl_sam-write_bam \
            --pe_frag_dist_std_dev 20 \
            --pe_frag_dist_mean 200 \
            --i-type fasta
        cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
        rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
    done
done

err_file="/dev/null"

for lc in se pe mp; do
    for mode in trans template; do
        for parser in memory stream; do
            # No, do not use SAM output.
            for i_type in fasta pbsim3_transcripts; do
                "${ART}" \
                    --builtin_qual_file HiSeq2500_125bp \
                    --i-file "${err_file}" \
                    --read_len 125 \
                    --mode "${mode}" \
                    --lc "${lc}" \
                    --i-parser "${parser}" \
                    --i-fcov "${FCOV}" \
                    --parallel "${PARALLEL}" \
                    --ins_rate_1 "${IDRATE}" \
                    --del_rate_1 "${IDRATE}" \
                    --o-fastq "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq \
                    --o-fasta "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fasta \
                    --o-pwa "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.pwa \
                    --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.hl.bam \
                    --o-hl_sam-num_threads 2 \
                    --o-hl_sam-compress_level u \
                    --o-hl_sam-write_bam \
                    --pe_frag_dist_std_dev 20 \
                    --pe_frag_dist_mean 200 \
                    --i-type "${i_type}"
                cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
                rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
            done
        done
    done
done

parser=memory
i_type=fasta
for lc in se pe mp; do
    for mode in wgs trans template; do
        # No, do not use SAM output.
        "${ART}" \
            --builtin_qual_file HiSeq2500_125bp \
            --i-file "${err_file}" \
            --read_len 125 \
            --mode "${mode}" \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov "${FCOV}" \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-fastq "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq \
            --o-fasta "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fasta \
            --o-pwa "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.pwa \
            --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.hl.bam \
            --o-hl_sam-num_threads 2 \
            --o-hl_sam-compress_level u \
            --o-hl_sam-write_bam \
            --pe_frag_dist_std_dev 20 \
            --pe_frag_dist_mean 200 \
            --i-type "${i_type}"
        cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
        rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
    done
done

unset FCOV
