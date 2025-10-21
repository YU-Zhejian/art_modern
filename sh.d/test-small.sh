#!/usr/bin/env bash
# shellcheck disable=SC1091
set -ue
SHDIR="$(readlink -f "$(dirname "${0}")")"
cd "${SHDIR}/../"

if [ ! -f data/raw_data/ce11.mRNA_head.cov_stranded.tsv ]; then
    python sh.d/test-small.sh.d/gen_cov.py
fi

export MPI_PARALLEL="${MPI_PARALLEL:-4}"
export SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-16}"
export IDRATE=0.1 # Increase indel rate to fail faster
OUT_DIR="$(readlink -f "opt/tmp/")"
export OUT_DIR
export ART="${ART:-opt/build_debug-install/art_modern}"
export MRNA_HEAD="data/raw_data/ce11.mRNA_head.fa"
export MRNA_PBSIM3_TRANSCRIPT="data/raw_data/ce11.mRNA_head.pbsim3.transcript"
export LAMBDA_PHAGE="data/raw_data/lambda_phage.fa"
export CE11_CHR1="data/raw_data/ce11_chr1.fa"
if [ -z "${MPIRUN:-}" ]; then
    ART_CMD_ASSEMBLED=("${ART}")
    export PARALLEL="4" # Reduce parallelism overhead for small tests
else
    export MPIRUN
    export PARALLEL="2"
    ART_CMD_ASSEMBLED=("${MPIRUN}" -np "${MPI_PARALLEL}" "${ART}")
fi

echo "ART=${ART} MPIRUN=${MPIRUN}"

function sam2bam() {
    # Single-threaded sorting should be fast enough
    samtools sort --write-index "${1}".sam -o "${1}".bam
    python sh.d/test-small.sh.d/test_sam.py "${2}" "${1}".bam
    rm -f "${1}".sam "${1}".bam "${1}".bam.csi "${1}".bam.bai
}

# TODO: Only used in 0_out_fmts.sh for checking.
function merge_file() {
    # Given $1: "${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq
    # Find:
    #   (With MPI)    Files like "${OUT_DIR}"/test_small_se_wgs_memory_sep.*.fastq
    #   (Without MPI) Files like "${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq
    # Merge into: "${OUT_DIR}"/test_small_se_wgs_memory_sep.fastq
    if [ -f "${1}" ]; then
        # No MPI run, nothing to do
        return
    fi
    # Only FASTQ, FASTA or SAM files can be merged. Get the extension.
    ext="${1##*.}"
    base="${1%.*}"
    files_to_merge=("${base}".*."${ext}")
    if [ "${ext}" == "fastq" ] || [ "${ext}" == "fq" ] || [ "${ext}" == "fasta" ] || [ "${ext}" == "fa" ]; then
        cat "${files_to_merge[@]}" >"${1}"
    elif [ "${ext}" == "sam" ] || [ "${ext}" == "bam" ]; then
        # Sort all files before merging
        new_files_to_merge=()
        for fn in "${files_to_merge[@]}"; do
            samtools sort --threads "${SAMTOOLS_THREADS}" --write-index "${fn}" -o "${fn}".sorted.bam
            new_files_to_merge+=("${fn}".sorted.bam)
        done
        samtools merge --threads "${SAMTOOLS_THREADS}" -o "${1}".sorted.bam "${new_files_to_merge[@]}"
        # Convert back to SAM if needed
        samtools view -h -o "${1}" "${1}".sorted.bam
        rm -fr "${1}".sorted.bam "${1}".sorted.bam.csi "${1}".sorted.bam.bai
        # Clean up sorted BAM files
        for fn in "${new_files_to_merge[@]}"; do
            rm -fr "${fn}" "${fn}".csi "${fn}".bai
        done
    else
        echo "merge_file: Unsupported extension: ${ext}" >&2
        exit 1
    fi
    # Remove original files
    for fn in "${files_to_merge[@]}"; do
        rm -f "${fn}"
    done
}

# Ensure OUT_DIR is clean
function assert_cleandir() {
    rm -d "${OUT_DIR}"
    mkdir "${OUT_DIR}"
}

function AM_EXEC() {
    echo "EXEC: ${ART_CMD_ASSEMBLED[*]} $*"
    "${ART_CMD_ASSEMBLED[@]}" "$@" &>>"${OUT_DIR}"/am_exec.log
    if [ ${?} -ne 0 ]; then
        echo "AM_EXEC failed with exit code ${?}" >&2
        cat "${OUT_DIR}"/am_exec.log >&2
        exit 1
    else
        echo "AM_EXEC succeeded."
        rm -f "${OUT_DIR}"/am_exec.log
    fi
    return ${?}
}

rm -fr "${OUT_DIR}" # Remove previous runs
mkdir "${OUT_DIR}"
. sh.d/test-small.sh.d/0_out_fmts.sh       # Test all output is working
. sh.d/test-small.sh.d/1_fail.sh           # FASTA that would fail the simulator
. sh.d/test-small.sh.d/2_wgs.sh            # WGS mode (with constant coverage)
. sh.d/test-small.sh.d/3_trans_constcov.sh # Transcript mode with constant coverage
. sh.d/test-small.sh.d/4_tmpl_constcov.sh  # Template mode with constant coverage
. sh.d/test-small.sh.d/5_trans_scov.sh     # Transcript mode with stranded/strandless coverage
. sh.d/test-small.sh.d/6_tmpl_scov.sh      # Template mode with stranded/strandless coverage
. sh.d/test-small.sh.d/7_tmpl_pbsim3.sh    # Transcript mode with pbsim3-formatted coverage
. sh.d/test-small.sh.d/8_trans_pbsim3.sh   # Template mode with pbsim3-formatted coverage
rm -d "${OUT_DIR}"                         # Which should now be empty
