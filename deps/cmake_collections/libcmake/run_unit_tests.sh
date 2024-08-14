#!/usr/bin/env bash
# Shell script that executes unit test.
# TODO: docs
mkdir -p build log
for fn in *_unittest; do
    printf "%s" "Testing ${fn}..."
    cd build &&
        rm -rf ./* &&
        cmake ../"${fn}"/ --trace --warn-uninitialized -Werror --log-level=DEBUG 2>&1 |
        tee ../log/"${fn}".log |
            grep '^--' >../log/"${fn}".status.log &&
        cd .. &&
        echo "SUCCESS" || echo "FAILED"
done
