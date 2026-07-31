#!/usr/bin/env bash
# Shell script that formats every source file tracked by Git.
set -ue

if [ "${PIXI_ENVIRONMENT_NAME:-}" = "fmt" ]; then
    true
else
    echo "This script should only be run in the 'fmt' environment."
    exit 1
fi

# Checked up front: a missing formatter must not look like a clean run.
MISSING=0
for tool in shfmt black clang-format cmake-format dos2unix; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "Required formatter not found: ${tool}" >&2
        MISSING=1
    fi
done
if [ "${MISSING}" -ne 0 ]; then
    exit 1
fi

NPROC="$(nproc)"

# Runs "${@:2}" over the file list on stdin, NPROC at a time; tracked-but-deleted files are skipped.
# Not '&' inside a 'while read' loop: that loop is a pipeline subshell, so a following 'wait' has no
# children and the stages overlap. xargs also fails loudly, where a bare 'wait' returns 0.
function run_formatter() {
    local label="${1}"
    shift
    local listed=()
    local files=()
    local f
    mapfile -t listed
    for f in "${listed[@]}"; do
        if [ -e "${f}" ]; then
            files+=("${f}")
        fi
    done
    if [ "${#files[@]}" -eq 0 ]; then
        return 0
    fi
    for f in "${files[@]}"; do
        echo "${label} ${f}"
    done
    printf '%s\0' "${files[@]}" | xargs -0 -P "${NPROC}" -n 8 "${@}"
}

git ls-files |
    grep -v 'deps' |
    grep -v '.idea/' |
    grep '\.sh$' |
    run_formatter SHFMT shfmt -i 4 -w

git ls-files |
    grep -v 'deps' |
    grep -v '.idea/' |
    grep -e '\.py$' -e '\.pyi$' |
    run_formatter BLACK black --line-length 120 --target-version py39 --quiet

git ls-files |
    grep -v 'deps/concurrentqueue/' |
    grep -v 'deps/labw_slim_htslib/' |
    grep -v 'deps/slim_fmt/' |
    grep -v 'deps/thread-pool/' |
    grep -v '.idea/' |
    grep -v 'benchmark_other_simulators/src/' |
    grep -e '\.cc$' -e '\.cpp$' -e '\.c$' -e '\.h$' -e '\.hh' -e '\.hpp$' |
    run_formatter CLANG-FORMAT clang-format -i

git ls-files |
    grep -v '.idea/' |
    grep -e '\.cmake$' -e 'CMakeLists.txt' |
    run_formatter CMAKE-FORMAT cmake-format --in-place

# Must stay last: it touches every tracked file, including those the formatters above rewrite.
git ls-files |
    run_formatter DOS2UNIX dos2unix --keepdate -q
