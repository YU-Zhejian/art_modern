#!/usr/bin/env bash
set -uex
cd "$(readlink -f "$(dirname "${0}")")"

if [ ! -f raw_data/ce11.mRNA_head.cov_stranded.tsv ]; then
    python test_small.sh.d/gen_cov.py
fi

. test_small.sh.d/base.sh
. test_small.sh.d/out_fmts.sh       # Test all output is working
. test_small.sh.d/wgs.sh            # WGS mode (with constant coverage)
. test_small.sh.d/trans_constcov.sh # Transcript mode with constant coverage
. test_small.sh.d/tmpl_constcov.sh  # Template mode with constant coverage
. test_small.sh.d/trans_scov.sh     # Transcript mode with stranded/strandless coverage
. test_small.sh.d/tmpl_scov.sh      # Template mode with stranded/strandless coverage
. test_small.sh.d/trans_pbsim3.sh   # Transcript mode with pbsim3-formatted coverage
. test_small.sh.d/tmpl_pbsim3.sh    # Template mode with pbsim3-formatted coverage
