#!/usr/bin/env bash
# Shell script that runs static analysis over Python and Shell sources.
set -u

RETV=0

# Compilation database used by clang-tidy. Shared with clangd; see .clangd.
BUILD_DIR="${AM_LINT_BUILD_DIR:-cmake-build-clangd}"

if which run-clang-tidy &>/dev/null && [ -e "${BUILD_DIR}/compile_commands.json" ]; then
    echo CLANG-TIDY
    # -Qunused-arguments: the database is generated for GCC and carries flags clang does not consume.
    run-clang-tidy -p "${BUILD_DIR}" -quiet -j "$(nproc)" \
        -extra-arg=-Qunused-arguments "^$(pwd)/src/.*" || RETV=1
fi

if which pyright &>/dev/null; then
    echo PYRIGHT
    pyright || RETV=1
fi

if which shellcheck &>/dev/null; then
    mapfile -t FILES < <(git ls-files | grep -v 'deps' | grep -v '.idea/' | grep '\.sh$')
    for line in "${FILES[@]}"; do
        if [ -e "${line}" ]; then
            echo SHELLCHECK "${line}"
            shellcheck "${line}" || RETV=1
        fi
    done
fi

exit "${RETV}"
