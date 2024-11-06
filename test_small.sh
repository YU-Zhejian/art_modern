#!/usr/bin/env bash
set -uex
cd "$(readlink -f "$(dirname "${0}")")"

if [ ! -f raw_data/ce11.mRNA_head.cov_stranded.tsv ]; then
    python test_small.sh.d/gen_cov.py
fi

. test_small.sh.d/base.sh
. test_small.sh.d/wgs.sh            # WGS mode (with constant coverage)
. test_small.sh.d/trans_constcov.sh # Transcript mode with constant coverage
. test_small.sh.d/tmpl_constcov.sh  # Template mode with constant coverage

# TODO: Trans mode with strandless coverage

# TODO: Trans mode with stranded coverage

# TODO: Trans mode with pbsim3-formatted coverage

# TODO: Template mode with strandless coverage

# TODO: Template mode with stranded coverage

# TODO: Template mode with pbsim3-formatted coverage
