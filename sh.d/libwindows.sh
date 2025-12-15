#shellcheck shell=bash

function is_windows() {
    UNAME_S_OUT="$(uname -s)"
    if [[ "${UNAME_S_OUT}" == *"_NT"* || "${UNAME_S_OUT}" == "MINGW"* || "${UNAME_S_OUT}" == "CYGWIN"* || "${UNAME_S_OUT}" == "MSYS"* ]]; then
        return 0
    else
        return 1
    fi
}
