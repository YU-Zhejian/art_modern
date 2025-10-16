#!/usr/bin/env bash
# shellcheck disable=SC1091
set -uex
cd "$(readlink -f "$(dirname "${0}")/../")"

if [ ! -f data/raw_data/ce11.mRNA_head.cov_stranded.tsv ]; then
    python sh.d/test-small.sh.d/gen_cov.py
fi

export PARALLEL="5" # Reduce parallelism overhead for small tests
export IDRATE=0.1   # Increase indel rate to fail faster
export OUT_DIR="$(readlink -f "opt/tmp/")"
export ART="${ART:-opt/build_debug/art_modern}"
export MRNA_HEAD="data/raw_data/ce11.mRNA_head.fa"
export MRNA_PBSIM3_TRANSCRIPT="data/raw_data/ce11.mRNA_head.pbsim3.transcript"
export LAMBDA_PHAGE="data/raw_data/lambda_phage.fa"

echo "ART=${ART}"

function sam2bam() {
    # Single-threaded sorting should be fast enough
    samtools sort --write-index "${1}".sam -o "${1}".bam
    python sh.d/test-small.sh.d/test_sam.py "${2}" "${1}".bam
    rm -f "${1}".sam "${1}".bam "${1}".bam.csi "${1}".bam.bai
}

rm -fr "${OUT_DIR}" # Remove previous runs
# . sh.d/test-small.sh.d/out_fmts.sh       # Test all output is working
# . sh.d/test-small.sh.d/fail.sh           # FASTA that would fail the simulator
. sh.d/test-small.sh.d/wgs.sh            # WGS mode (with constant coverage)
. sh.d/test-small.sh.d/trans_constcov.sh # Transcript mode with constant coverage
. sh.d/test-small.sh.d/tmpl_constcov.sh  # Template mode with constant coverage
. sh.d/test-small.sh.d/trans_scov.sh     # Transcript mode with stranded/strandless coverage
. sh.d/test-small.sh.d/tmpl_scov.sh      # Template mode with stranded/strandless coverage
. sh.d/test-small.sh.d/trans_pbsim3.sh   # Transcript mode with pbsim3-formatted coverage
. sh.d/test-small.sh.d/tmpl_pbsim3.sh    # Template mode with pbsim3-formatted coverage
rm -d "${OUT_DIR}"                       # Which should now be empty
