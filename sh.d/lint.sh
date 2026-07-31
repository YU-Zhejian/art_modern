#!/usr/bin/env bash
# Shell script that runs static analysis over Python and Shell sources.
set -u

RETV=0

# Compilation database used by clang-tidy. Shared with clangd; see .clangd.
BUILD_DIR="${AM_LINT_BUILD_DIR:-cmake-build-clangd}"

# Checked up front: a missing tool must not look like a lint-clean tree.
MISSING=0
for tool in run-clang-tidy pyright shellcheck; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "Required tool not found: ${tool}" >&2
        MISSING=1
    fi
done
if [ ! -e "${BUILD_DIR}/compile_commands.json" ]; then
    echo "Compilation database not found: ${BUILD_DIR}/compile_commands.json" >&2
    echo "Run 'make lint', which configures it, or point AM_LINT_BUILD_DIR at a build directory" \
        "configured with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON." >&2
    MISSING=1
fi
if [ "${MISSING}" -ne 0 ]; then
    exit 1
fi

echo CLANG-TIDY
# -Qunused-arguments: the database is generated for GCC and carries flags clang does not consume.
run-clang-tidy -p "${BUILD_DIR}" -quiet -j "$(nproc)" \
    -extra-arg=-Qunused-arguments "^$(pwd)/src/.*" || RETV=1

echo PYRIGHT
pyright || RETV=1

mapfile -t FILES < <(git ls-files | grep -v 'deps' | grep -v '.idea/' | grep '\.sh$')
for line in "${FILES[@]}"; do
    if [ -e "${line}" ]; then
        echo SHELLCHECK "${line}"
        shellcheck "${line}" || RETV=1
    fi
done

exit "${RETV}"
