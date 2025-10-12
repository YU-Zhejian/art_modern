#!/bin/bash
set -ue
mkdir -p opt
# MKL is not tested.

BUILD_ONLY_TEST=1
USE_RANDOM_GENERATOR=STL
USE_QUAL_GEN=STL
for CMAKE_BUILD_TYPE in Debug Release RelWithDebInfo; do
    for USE_THREAD_PARALLEL in ASIO BS NOP; do
        for USE_MALLOC in NOP JEMALLOC MIMALLOC; do
            for USE_LIBFMT in "UNSET" "fmt"; do
                for USE_HTSLIB in "UNSET" "hts"; do
                    for USE_CONCURRENT_QUEUE in "UNSET" "SYS"; do
                        for USE_ABSL in "UNSET" "SYS"; do
                            CMAKE_FLAGS=""
                            CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                            CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_RANDOM_GENERATOR=${USE_RANDOM_GENERATOR}"
                            CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_THREAD_PARALLEL=${USE_THREAD_PARALLEL}"
                            CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_QUAL_GEN=${USE_QUAL_GEN}"
                            CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_MALLOC=${USE_MALLOC}"
                            if [ "${USE_LIBFMT}" != "UNSET" ]; then
                                CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_LIBFMT=${USE_LIBFMT}"
                            fi
                            if [ "${USE_HTSLIB}" != "UNSET" ]; then
                                CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_HTSLIB=${USE_HTSLIB}"
                            fi
                            if [ "${USE_ABSL}" != "UNSET" ]; then
                                CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_ABSL=SYS"
                            fi
                            if [ "${USE_CONCURRENT_QUEUE}" != "UNSET" ]; then
                                if [ -f /usr/include/concurrentqueue/moodycamel/concurrentqueue.h ]; then
                                    CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_CONCURRENT_QUEUE=/usr/include/concurrentqueue/moodycamel"
                                elif [ -f /usr/include/concurrentqueue/concurrentqueue.h ]; then
                                    CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_CONCURRENT_QUEUE=/usr/include/concurrentqueue"
                                else
                                    echo "concurrentqueue.h not found!"
                                    exit 1
                                fi
                            fi
                            echo "Testing with CMAKE_FLAGS: ${CMAKE_FLAGS}" | tee opt/testbuild.log
                            make testbuild-child \
                                CMAKE_FLAGS="${CMAKE_FLAGS}" \
                                BUILD_ONLY_TEST="${BUILD_ONLY_TEST}" \
                                &>>opt/testbuild.log || exit 1
                        done
                    done
                done
            done
        done
    done
done
BUILD_ONLY_TEST=2
CMAKE_BUILD_TYPE=RelWithDebInfo
for USE_RANDOM_GENERATOR in STL PCG GSL BOOST; do
    for USE_QUAL_GEN in WALKER STL; do
        CMAKE_FLAGS=""
        CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_RANDOM_GENERATOR=${USE_RANDOM_GENERATOR}"
        CMAKE_FLAGS="${CMAKE_FLAGS} -DUSE_QUAL_GEN=${USE_QUAL_GEN}"

        echo "Testing with CMAKE_FLAGS: ${CMAKE_FLAGS}" | tee opt/testbuild.log
        make testbuild-child \
            CMAKE_FLAGS="${CMAKE_FLAGS}" \
            BUILD_ONLY_TEST="${BUILD_ONLY_TEST}" \
            &>>opt/testbuild.log || exit 1
    done
done
