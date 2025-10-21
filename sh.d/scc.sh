#!/usr/bin/env bash
# Shell script that counts how many lines of code was written by us.
# shellcheck disable=SC2086

builtin set -ue
NAME="scc.sh"
VERSION=0.1
SCC="scc"

for requested_binaries in git jq scc awk grep sed xargs; do
    if builtin type -p "${requested_binaries}" &>/dev/null; then
        true
    else
        echo "${requested_binaries} required!"
        builtin exit 1
    fi
done

echo "Enumerating sources ..."

SOURCES=$(
    git ls-files |
        grep -v '\.idea' |
        grep -v 'deps' |
        grep -v 'Illumina_profiles' |
        grep -v 'benchmark_other_simulators/src' |
        xargs
)
echo "Enumerating sources FIN"

CLOC_INFO=$("${SCC}" ${SOURCES})

cat <<EOF
${NAME} ver. ${VERSION}
Called by: ${0} ${*}
Code count:
${CLOC_INFO}
EOF

# Inflation Detection
echo "Detecting external C/C++ code ..."

for dir in deps/* explore/*; do
    if [ ! -d "${dir}" ]; then
        continue
    fi
    SOURCES=$(git ls-files "${dir}" | xargs)
    "${SCC}" --include-ext c,cc,h,hh,cpp,hpp,cxx,hxx --format json ${SOURCES} >deps_code_count.json
    COUNTS=$(jq 'reduce .[] as $item (0; . + $item.Count)' deps_code_count.json)
    LINES=$(jq 'reduce .[] as $item (0; . + $item.Lines)' deps_code_count.json)
    printf "Directory: %-48s has %-10d files with %-6d lines of C/C++ code\n" "${dir}" "${COUNTS}" "${LINES}"
    rm -f deps_code_count.json
done
builtin exit 0
