#shellcheck shell=bash

function is_windows(){
    if [[ "$(uname -s)" == *"_NT"* || "$(uname -s)" == "MINGW"* || "$(uname -s)" == "CYGWIN"* || "$(uname -s)" == "MSYS"* ]]; then
        return 0
    else
        return 1
    fi
}
