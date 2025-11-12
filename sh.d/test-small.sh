#!/usr/bin/env bash
# shellcheck disable=SC1091
# Taking environment variables:
# - NO_FASTQC: If set to "1", skip fastqc tests
# - HELP_VERSION_ONLY: If set to "1", only test help and version outputs
# - FORMAT_ONLY: If set to "1", only test output formats
# - MPIEXEC: If set, use this command to run ART in MPI mode
# - ART: If set, use this path as the art_modern executable
# - MPI_PARALLEL: If set, use this many MPI processes
# - SAMTOOLS_THREADS: If set, use this many threads for samtools operations

set -ue
SHDIR="$(readlink -f "$(dirname "${0}")")"
cd "${SHDIR}/../"

if [ ! -f data/raw_data/ce11.mRNA_head.cov_stranded.tsv ]; then
    python "${SHDIR}"/test-small.sh.d/gen_cov.py data/raw_data/ce11.mRNA_head 5
fi

# Create a temporary output directory
OUT_DIR="$(mktemp -d art_modern_test_small.d.XXXXXX --tmpdir=/tmp)"

MRNA_HEAD="$(readlink -f "data/raw_data/ce11.mRNA_head.fa")"
MRNA_PBSIM3_TRANSCRIPT="$(readlink -f "data/raw_data/ce11.mRNA_head.pbsim3.transcript")"
LAMBDA_PHAGE="$(readlink -f "data/raw_data/lambda_phage.fa")"
CE11_CHR1="$(readlink -f "data/raw_data/ce11_chr1.fa")"

export MPI_PARALLEL="${MPI_PARALLEL:-4}"
export SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-16}" # Also used by fastqc
export IDRATE=0.1                                 # Increase indel rate to fail faster
export ART_MODERN_PATH="${ART_MODERN_PATH}"       # Do NOT have a default, must be set from outside
export APB_PATH="${APB_PATH}"                     # Do NOT have a default, must be set from outside
export LAMBDA_PHAGE
export CE11_CHR1
export OUT_DIR
export MRNA_HEAD
export MRNA_PBSIM3_TRANSCRIPT

# Assemble the ART command
# If MPIEXEC is set, use it to run ART in parallel using MPI
if [ -z "${MPIEXEC:-}" ]; then
    ART_CMD_ASSEMBLED=("${ART_MODERN_PATH}")
    APB_CMD_ASSEMBLED=("${APB_PATH}")
    export PARALLEL="${PARALLEL:-4}" # Reduce parallelism overhead for small tests
else
    export MPIEXEC
    export PARALLEL="${PARALLEL:-2}"
    ART_CMD_ASSEMBLED=("${MPIEXEC}" -n "${MPI_PARALLEL}" "${ART_MODERN_PATH}")
    APB_CMD_ASSEMBLED=("${MPIEXEC}" -n "${MPI_PARALLEL}" "${APB_PATH}")
fi

echo "ART_MODERN_PATH=${ART_MODERN_PATH} MPIEXEC=${MPIEXEC} OUT_DIR=${OUT_DIR}"

function sam2bam() {
    # Single-threaded sorting should be fast enough
    samtools sort --write-index "${1}".sam -o "${1}".bam
    python sh.d/test-small.sh.d/test_sam.py "${2}" "${1}".bam
    rm -f "${1}".sam "${1}".bam "${1}".bam.csi "${1}".bam.bai
}

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
    rm -fr "${OUT_DIR}"/log_*.d "${OUT_DIR}"/*.log # Remove logs
    rm -d "${OUT_DIR}"
    mkdir "${OUT_DIR}"
}

EXEC_ORDER=0

function AM_EXEC() {
    echo "EXEC ${EXEC_ORDER}: ${ART_CMD_ASSEMBLED[*]} $*"
    env \
        "ART_LOG_DIR=${OUT_DIR}/log_${EXEC_ORDER}.d" \
        "${ART_CMD_ASSEMBLED[@]}" "$@" &>>"${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log
    retval=${?}
    if [ ${retval} -ne 0 ]; then
        echo "AM_EXEC failed with exit code ${retval}" >&2
        cat "${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log >&2
        exit 1
    else
        echo "AM_EXEC succeeded."
        rm -f "${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log
    fi
    EXEC_ORDER=$((EXEC_ORDER + 1))
    return ${retval}
}
function APB_EXEC() {
    echo "EXEC ${EXEC_ORDER}: ${APB_CMD_ASSEMBLED[*]} $*"
    env \
        "ART_LOG_DIR=${OUT_DIR}/log_${EXEC_ORDER}.d" \
        "${APB_CMD_ASSEMBLED[@]}" "$@" &>>"${OUT_DIR}"/apb_exec_"${EXEC_ORDER}".log
    retval=${?}
    if [ ${retval} -ne 0 ]; then
        echo "APB_EXEC failed with exit code ${retval}" >&2
        cat "${OUT_DIR}"/apb_exec_"${EXEC_ORDER}".log >&2
        exit 1
    else
        echo "APB_EXEC succeeded."
        rm -f "${OUT_DIR}"/apb_exec_"${EXEC_ORDER}".log
    fi
    EXEC_ORDER=$((EXEC_ORDER + 1))
    return ${retval}
}

rm -fr "${OUT_DIR}" # Remove previous runs
mkdir "${OUT_DIR}"
AM_EXEC --version  # Just to log the version used
AM_EXEC --help     # Just to log the help message
APB_EXEC --version # Just to log the version used
APB_EXEC --help    # Just to log the help message
if [ "${HELP_VERSION_ONLY:-0}" == "1" ]; then
    rm -fr "${OUT_DIR}"
    exit 0
fi
. "${SHDIR}"/test-small.sh.d/0_out_fmts.sh       # Test all output is working
. "${SHDIR}"/test-small.sh.d/1_fail.sh           # FASTA that would fail the simulator
. "${SHDIR}"/test-small.sh.d/2_wgs.sh            # WGS mode (with constant coverage)
. "${SHDIR}"/test-small.sh.d/3_trans_constcov.sh # Transcript mode with constant coverage
. "${SHDIR}"/test-small.sh.d/4_tmpl_constcov.sh  # Template mode with constant coverage
. "${SHDIR}"/test-small.sh.d/5_trans_scov.sh     # Transcript mode with stranded/strandless coverage
. "${SHDIR}"/test-small.sh.d/6_tmpl_scov.sh      # Template mode with stranded/strandless coverage
. "${SHDIR}"/test-small.sh.d/7_tmpl_pbsim3.sh    # Transcript mode with pbsim3-formatted coverage
. "${SHDIR}"/test-small.sh.d/8_trans_pbsim3.sh   # Template mode with pbsim3-formatted coverage
. "${SHDIR}"/test-small.sh.d/21-apb-se.sh        # APB single-end test
. "${SHDIR}"/test-small.sh.d/22-apb-pe.sh        # APB paired-end test
rm -d "${OUT_DIR}"                               # Which should now be empty
