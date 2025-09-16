#!/bin/bash
set -ue
mkdir -p opt
# MKL is not tested.
for CMAKE_BUILD_TYPE in Debug Release RelWithDebInfo SomeNonesense; do
    for USE_RANDOM_GENERATOR in STL PCG GSL BOOST; do
        for USE_THREAD_PARALLEL in ASIO BS NOP; do
            for USE_QUAL_GEN in WALKER STL; do
                for USE_MALLOC in NOP JEMALLOC MIMALLOC; do
                    CMAKE_FLAGS=""
                    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                    CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_RANDOM_GENERATOR=${USE_RANDOM_GENERATOR}"
                    CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_THREAD_PARALLEL=${USE_THREAD_PARALLEL}"
                    CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_QUAL_GEN=${USE_QUAL_GEN}"
                    CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_MALLOC=${USE_MALLOC}"

                    echo "Testing with CMAKE_FLAGS: ${CMAKE_FLAGS}" | tee opt/testbuild.log
                    make testbuild-child CMAKE_FLAGS="${CMAKE_FLAGS}" &>> opt/testbuild.log || exit 1
                done
            done
        done
    done
done
