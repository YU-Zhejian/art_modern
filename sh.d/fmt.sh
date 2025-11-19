#!/usr/bin/env bash
# Shell script that counts how many lines of code was written by us.
# shellcheck disable=SC2086

which shfmt &>/dev/null &&
    git ls-files |
    grep -v 'deps' |
        grep -v '.idea/' |
        grep '\.sh$' |
        while read -r line; do
            if [ -e "${line}" ]; then
                {
                    echo SHFMT "${line}"
                    shfmt -i 4 -w "${line}"
                } &
            fi
        done
wait

which black &>/dev/null &&
    git ls-files |
    grep -v 'deps' |
        grep -v '.idea/' |
        grep -e '\.py$' -e '\.pyi$' |
        while read -r line; do
            if [ -e "${line}" ]; then
                {
                    echo BLACK "${line}"
                    black --line-length 120 --target-version py39 --quiet "${line}"
                } &
            fi
        done
wait

which clang-format &>/dev/null &&
    git ls-files |
    grep -v 'deps/concurrentqueue/' |
        grep -v 'deps/labw_slim_htslib/' |
        grep -v 'deps/slim_fmt/' |
        grep -v 'deps/thread-pool/' |
        grep -v '.idea/' |
        grep -v 'benchmark_other_simulators/src/' |
        grep -e '\.cc$' -e '\.cpp$' -e '\.c$' -e '\.h$' -e '\.hh' -e '\.hpp$' |
        while read -r line; do
            if [ -e "${line}" ]; then
                {
                    echo "CLANG-FORMAT ${line}"
                    clang-format -i "${line}"
                } &
            fi
        done
wait

which cmake-format &>/dev/null &&
    git ls-files |
    grep -v '.idea/' |
        grep -e '\.cmake$' -e 'CMakeLists.txt' |
        while read -r line; do
            if [ -e "${line}" ]; then
                {
                    echo "CMAKE-FORMAT ${line}"
                    cmake-format \
                        --tab-size 4 \
                        --line-width 120 \
                        --in-place "${line}"
                } &
            fi
        done
wait

which dos2unix &>/dev/null &&
    git ls-files |
    while read -r line; do
        if [ -e "${line}" ]; then
            {
                echo DOS2UNIX "${line}"
                dos2unix --keepdate "${line}" &>>/dev/null
            } &
        fi
    done
wait
