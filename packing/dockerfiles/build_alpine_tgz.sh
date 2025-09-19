#!/usr/bin/env sh
set -ue
# Get mounted paths
DIST_TGZ_PATH="/mnt/art_modern-master"

# Set locales
export LANG=C
export LC_ALL=C
export LC_CTYPE=C
export LC_NUMERIC=C
export LC_TIME=C
export LC_COLLATE=C
export LC_MONETARY=C
export LC_MESSAGES=C
export LC_PAPER=C
export LC_NAME=C
export LC_ADDRESS=C
export LC_TELEPHONE=C
export LC_MEASUREMENT=C
export LC_IDENTIFICATION=C

cd "${DIST_TGZ_PATH}"
make rel_with_dbg_alpine
