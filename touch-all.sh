#!/usr/bin/env bash
find . | grep -v 'deps' | while read -r fn; do
    echo TOUCH "${fn}"
    touch "${fn}" &
done
wait
