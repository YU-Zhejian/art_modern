#!/usr/bin/env bash
set -ue
git ls-files |
    grep -v '\.idea' |
    grep -v 'deps' |
    grep -v 'Illumina_profiles' |
    grep -v 'benchmark_other_simulators/src' |
    grep -e '\.cc$' -e '\.cpp$' -e '\.cxx$' -e '\.c$' -e '\.h$' -e '\.hpp$' -e '\.hxx$' -e '\.hh$' |
    while read -r f; do
        cat "${f}" |
            sed -e 's;^\s+;;' |
            grep '^#include' |
            grep -e '<boost/' -e '"boost/' |
            sed -e 's;[[:space:]]*//.*;;' |
            sed 's;$; // '"${f}"';'
    done |
    sort |
    uniq
