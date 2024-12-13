#!/usr/bin/env bash
set -uex
cd "$(readlink -f "$(dirname "${0}")/../")"

if [ ! -f data/raw_data/ce11.mRNA_head.cov_stranded.tsv ]; then
    python sh.d/test_small.sh.d/gen_cov.py
fi

export PARALLEL=0
export IDRATE=0.1
export OUT_DIR="opt/tmp/"
export ART="opt/build_debug/art_modern"
export MRNA_HEAD="data/raw_data/ce11.mRNA_head.fa"
export ERR_FA="data/raw_data/err.fa"

function sam2bam() {
    samtools sort -@20 --write-index "${1}".sam -o "${1}".bam
    python sh.d/test_small.sh.d/test_sam.py "${2}" "${1}".bam
    rm -f "${1}".sam "${1}".bam "${1}".bam.csi "${1}".bam.bai
}

rm -fr "${OUT_DIR}"                      # Remove previous runs
. sh.d/test_small.sh.d/out_fmts.sh       # Test all output is working
. sh.d/test_small.sh.d/fail.sh           # FASTA that would fail the simulator
. sh.d/test_small.sh.d/wgs.sh            # WGS mode (with constant coverage)
. sh.d/test_small.sh.d/trans_constcov.sh # Transcript mode with constant coverage
. sh.d/test_small.sh.d/tmpl_constcov.sh  # Template mode with constant coverage
. sh.d/test_small.sh.d/trans_scov.sh     # Transcript mode with stranded/strandless coverage
. sh.d/test_small.sh.d/tmpl_scov.sh      # Template mode with stranded/strandless coverage
. sh.d/test_small.sh.d/trans_pbsim3.sh   # Transcript mode with pbsim3-formatted coverage
. sh.d/test_small.sh.d/tmpl_pbsim3.sh    # Template mode with pbsim3-formatted coverage
rm -d "${OUT_DIR}"                       # Which should now be empty
