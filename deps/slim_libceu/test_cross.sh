#!/usr/bin/env bash
set -uex
rm -fr opt
LINE_ID=0
sed '1d' cross_config.csv | grep -v '^#' |
    while read -r line; do
        LINE_ID=$((LINE_ID + 1))
        # Split the line using ,
        gcc_triplet=$(echo "${line}" | cut -d, -f1)
        qemu_postfix=$(echo "${line}" | cut -d, -f2)
        gcc_options=$(echo "${line}" | cut -d, -f3)
        export CXX="${gcc_triplet}-g++"
        export CC="${gcc_triplet}-gcc"
        export CXXFLAGS="${gcc_options}"
        export CFLAGS="${gcc_options}"
        mkdir -p opt/"${LINE_ID}"
        cmake -B opt/"${LINE_ID}" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_CXX_COMPILER="${CXX}" \
            -DCMAKE_C_COMPILER="${CC}" \
            -DBUILD_SHARED_LIBS=OFF \
            .
        cmake --build opt/"${LINE_ID}"
        qemu-"${qemu_postfix}" opt/"${LINE_ID}"/ceu_check_exe >opt/"${LINE_ID}"/test_output.txt 2>&1
    done
