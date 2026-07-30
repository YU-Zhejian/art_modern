#!/usr/bin/env bash
# Shell script that formats every source file tracked by Git.
set -ue

if [ "${PIXI_ENVIRONMENT_NAME:-}" = "fmt" ]; then
    true
else
    echo "This script should only be run in the 'fmt' environment."
    exit 1
fi

# A missing formatter must not be mistaken for a clean run, so check up front rather than skipping
# silently. All of these come from '[feature.fmt.dependencies]' in pixi.toml.
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

# Runs "${@:2}" over the file list read from stdin, at most NPROC invocations at a time. Files that
# Git still tracks but that no longer exist are skipped.
#
# Do not go back to backgrounding the formatters with '&' from inside a 'while read' loop: such a
# loop is the last element of a pipeline and therefore runs in a subshell, so a 'wait' placed after
# the pipeline has no children of its own to wait for. The stages then overlap, and dos2unix can
# rewrite a file while clang-format is still working on it. xargs also reports a non-zero status
# when an invocation fails, which a bare 'wait' never does.
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
