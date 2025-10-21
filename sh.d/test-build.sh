#!/bin/bash
set -ue
mkdir -p opt

if [ ! -z "${CMAKE_FLAGS:-}" ]; then
    OLD_CMAKE_FLAGS="${CMAKE_FLAGS}"
    echo "Old CMAKE_FLAGS: ${OLD_CMAKE_FLAGS}"
else
    OLD_CMAKE_FLAGS=""
fi

function do_build() {
    CMAKE_FLAGS="${OLD_CMAKE_FLAGS}"
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
}

SHDIR="$(dirname "$(readlink -f "${0}")")"

# A list of possible random generators
RANDOM_GENERATORS=(STL PCG BOOST)

if cmake -P "${SHDIR}/test-build.d/test-gsl.cmake" &>>/dev/null; then
    RANDOM_GENERATORS+=(GSL)
fi
if cmake -P "${SHDIR}/test-build.d/test-mkl.cmake" &>>/dev/null; then
    RANDOM_GENERATORS+=(MKL)
fi

# Default options
USE_RANDOM_GENERATOR=PCG
USE_QUAL_GEN=WALKER
USE_MALLOC=AUTO
USE_LIBFMT=UNSET
USE_HTSLIB=UNSET
USE_CONCURRENT_QUEUE=UNSET
USE_ABSL=UNSET

for CMAKE_BUILD_TYPE in Debug Release RelWithDebInfo; do
    BUILD_ONLY_TEST=1

    for USE_THREAD_PARALLEL in ASIO BS NOP; do
        do_build
    done
    USE_THREAD_PARALLEL=ASIO

    for USE_MALLOC in NOP JEMALLOC MIMALLOC; do
        do_build
    done
    USE_MALLOC=AUTO

    for USE_LIBFMT in "UNSET" "fmt"; do
        do_build
    done
    USE_LIBFMT=UNSET

    for USE_HTSLIB in "UNSET" "hts"; do
        do_build
    done
    USE_HTSLIB=UNSET

    for USE_CONCURRENT_QUEUE in "UNSET" "SYS"; do
        do_build
    done
    USE_CONCURRENT_QUEUE=UNSET

    for USE_ABSL in "UNSET" "SYS"; do
        do_build
    done
    USE_ABSL=UNSET

    BUILD_ONLY_TEST=2 # Now we requires a testsmall to be run

    for USE_RANDOM_GENERATOR in "${RANDOM_GENERATORS[@]}"; do
        do_build
    done
    USE_RANDOM_GENERATOR=PCG

    # TODO: Migrate this thing to unit tests

    for USE_QUAL_GEN in WALKER STL; do
        do_build
    done
    USE_QUAL_GEN=WALKER
done
