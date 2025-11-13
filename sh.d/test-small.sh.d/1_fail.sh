# shellcheck shell=bash
FCOV=10
parser=memory
err_file="${PROJDIR}/data/raw_data/err.fa"
for lc in se pe mp; do
    for mode in wgs trans template; do
        AM_EXEC \
            --i-file "${err_file}" \
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
        merge_file "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq
        cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
        rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
        assert_cleandir
    done
done

parser=htslib
mode=wgs
for lc in se pe mp; do
    AM_EXEC \
        --i-file "${err_file}" \
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
    merge_file "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq
    cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
    rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
    assert_cleandir
done

parser=stream
for lc in se pe mp; do
    for mode in trans template; do
        AM_EXEC \
            --i-file "${err_file}" \
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
        merge_file "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq
        cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
        rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
        assert_cleandir
    done
done

err_file="${OUT_DIR}/dev_null"

for lc in se pe mp; do
    for mode in trans template; do
        for parser in memory stream; do
            # No, do not use SAM output.
            for i_type in fasta pbsim3_transcripts; do
                touch "${err_file}"
                AM_EXEC \
                    --i-file "${err_file}" \
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
                merge_file "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq
                cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
                rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
                cmp "${err_file}" /dev/null
                rm -f "${err_file}"
                assert_cleandir
            done
        done
    done
done

parser=memory
i_type=fasta
for lc in se pe mp; do
    for mode in wgs trans template; do
        # No, do not use SAM output.
        touch "${err_file}"
        AM_EXEC \
            --i-file "${err_file}" \
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
        merge_file "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq
        cmp "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.fastq /dev/null
        rm -f "${OUT_DIR}"/test_small_"${lc}"_"${mode}"_"${parser}"_willfail.*
        cmp "${err_file}" /dev/null
        rm -f "${err_file}"
        assert_cleandir
    done
done

unset FCOV
