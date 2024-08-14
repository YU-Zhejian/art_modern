#!/usr/bin/env sh
# This file detects C/C++ preprocesor macros.
# Synopsis: ${0} ${CC} ${LANG} ${OUT}
# Params:
#   - ${CC} The C/C++ compiler, should be absolute path to `gcc`, `clang`, etc.
#   - ${LANG} Language, should be `c` or `c++`.
#   - ${OUT} The output directory.
# Output:
# If there's no error, there would be a file at "${CEU_DCPM_OUT}/test_${CEU_DCPM_LANG}_cpp.h" with all macros inside.
set -ue
CEU_DCPM_CC=${1}
CEU_DCPM_LANG=${2}
CEU_DCPM_OUT=${3}

mkdir -p "${CEU_DCPM_OUT}"

perform_cmd() {
    eval "${1}" >"${CEU_DCPM_OUT}/test_${CEU_DCPM_LANG}_cpp.h" \
        2>"${CEU_DCPM_OUT}/test_${CEU_DCPM_LANG}_cpp.err"
}

for CMD in \
    "\"${CEU_DCPM_CC}\" -dM -E -x \"${CEU_DCPM_LANG}\" - < /dev/null" \
    "\"${CEU_DCPM_CC}\" -dM -E -x \"${CEU_DCPM_LANG}\" /dev/null" \
    "\"${CEU_DCPM_CC}\" -dM -E -x\"${CEU_DCPM_LANG}\" - < /dev/null" \
    "\"${CEU_DCPM_CC}\" -dM -E -x\"${CEU_DCPM_LANG}\" < /dev/null" \
    "\"${CEU_DCPM_CC}\" -dM -E - < /dev/null" \
    "\"${CEU_DCPM_CC}\" -dM -E - /dev/null"; do
    if perform_cmd "${CMD}"; then
        exit
    fi
done
