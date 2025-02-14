#!/bin/bash
set -uex
# MKL is not tested.
for CMAKE_BUILD_TYPE in Debug Release RelWithDebInfo SomeNonesense; do
    for USE_RANDOM_GENERATOR in STL PCG GSL BOOST; do
        for USE_THREAD_PARALLEL in ASIO BS NOP; do
            for USE_QUAL_GEN in WALKER STL; do
                for USE_MALLOC in NOP JEMALLOC MIMALLOC; do
                    make testbuild-child CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DUSE_RANDOM_GENERATOR=${USE_RANDOM_GENERATOR} -DUSE_THREAD_PARALLEL=${USE_THREAD_PARALLEL} -DUSE_QUAL_GEN=${USE_QUAL_GEN} -DUSE_MALLOC=${USE_MALLOC}"
                done
            done
        done
    done
done
