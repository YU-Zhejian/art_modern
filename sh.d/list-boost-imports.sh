#!/usr/bin/env bash
set -ue
git ls-files |
    grep -v '\.idea' |
    grep -v 'deps' |
    grep -v 'Illumina_profiles' |
    grep -v 'benchmark_other_simulators/src' |
    grep -e '\.cc$' -e '\.c$' -e '\.hh' -e '\.hh' |
    while read -r f; do
        cat "${f}" |
            grep '^#include' |
            grep '<boost/' |
            sed -e 's;[[:space:]]*//.*;;' |
            sed 's;$; // '"${f}"';'
    done |
    sort |
    uniq
