# TODO: Check whether the architecture is ARM
# If so, disable MKL and VMT19937.
# Check for: __ARM_ARCH __arm __arm__
# See also: https://gitlab.arm.com/libraries/openrng
# See also: https://developer.arm.com/documentation/101004/2404/General-information/Access-Arm-Performance-Libraries?lang=en#access-pl
# See also: https://developer.arm.com/Tools%20and%20Software/Arm%20Compiler%20for%20Linux#Downloads

cmake_minimum_required(VERSION 3.17)

if(NOT DEFINED BOOST_CONFIG_PROVIDED_BY_BOOST)
    set(BOOST_CONFIG_PROVIDED_BY_BOOST ON)
endif()
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
    if(BOOST_CONFIG_PROVIDED_BY_BOOST)
        cmake_policy(SET CMP0167 NEW)
    else()
        cmake_policy(SET CMP0167 OLD)
    endif()
endif()

project(bench_rand C CXX)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(LINK_LIBS "")

# Find Boost
find_package(Boost REQUIRED COMPONENTS random filesystem)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    set(LINK_LIBS Boost::filesystem ${LINK_LIBS})
endif()

# Find MKL
find_package(MKL)
if(MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIRS})
    set(LINK_LIBS MKL::MKL ${LINK_LIBS})
endif()

# Find GSL
find_package(GSL)
if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    set(LINK_LIBS GSL::gsl ${LINK_LIBS})
endif()

# Find Abseil
find_package(absl REQUIRED)
if(absl_FOUND)
    include_directories(${absl_INCLUDE_DIRS})
    set(LINK_LIBS absl::random_random absl::random_distributions absl::base ${LINK_LIBS})
endif()

# Find ARMPL (ARM Performance Libraries)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ARMPL armpl)
if(ARMPL_FOUND)
    include_directories(${ARMPL_INCLUDEDIR})
    set(LINK_LIBS ${ARMPL_LIBRARIES} ${LINK_LIBS})
    add_definitions(-DARMPL_FOUND)
endif()

# Find TestU01. See <https://packages.debian.org/trixie/libtestu01-0-dev>
find_library(TESTU01_LIB NAMES testu01)
find_library(TESTU01_MYLIB NAMES testu01mylib mylib)
find_library(TESTU01_PROBDIST NAMES testu01probdist probdist)

include_directories("${CMAKE_CURRENT_LIST_DIR}/../../deps/pcg-cpp-0.98/include/")
include_directories("${CMAKE_CURRENT_LIST_DIR}/../../deps/VMT19937-master/include/")
include_directories("${CMAKE_CURRENT_LIST_DIR}/../../deps/xoshiro/include/")
include_directories("${CMAKE_CURRENT_LIST_DIR}/../../deps/other_rngs/")
include_directories("${CMAKE_CURRENT_LIST_DIR}/include/")

configure_file("${CMAKE_CURRENT_LIST_DIR}/include/bench_rand_conf.hh.in"
               "${CMAKE_CURRENT_BINARY_DIR}/include/bench_rand_conf.hh")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include/")

set(COMPILE_OPTIONS
    -O3
    -march=native
    -mtune=native
    -Wall
    -Wextra
    -Wpedantic
    -Og
    -g)

add_executable(bench_common_task bench_common_task.cc)
target_link_libraries(bench_common_task PUBLIC ${LINK_LIBS})
target_compile_options(bench_common_task PUBLIC ${COMPILE_OPTIONS})

add_executable(bench_bits bench_bits.cc)
target_link_libraries(bench_bits PUBLIC ${LINK_LIBS})
target_compile_options(bench_bits PUBLIC ${COMPILE_OPTIONS})

add_executable(bench_dists bench_dists.cc)
target_link_libraries(bench_dists PUBLIC ${LINK_LIBS})
target_compile_options(bench_dists PUBLIC ${COMPILE_OPTIONS})

add_executable(bench_uniform_real bench_uniform_real.cc)
target_link_libraries(bench_uniform_real PUBLIC ${LINK_LIBS})
target_compile_options(bench_uniform_real PUBLIC ${COMPILE_OPTIONS})

if(TESTU01_LIB_FOUND
   AND TESTU01_MYLIB_FOUND
   AND TESTU01_PROBDIST_FOUND)
    add_executable(bench_testu01 bench_testu01.cc)
    target_link_libraries(bench_testu01 PUBLIC ${LINK_LIBS} ${TESTU01_LIB} ${TESTU01_MYLIB} ${TESTU01_PROBDIST})
    target_compile_options(bench_testu01 PUBLIC ${COMPILE_OPTIONS})

    add_executable(testu01_main testu01_main.c)
    target_link_libraries(testu01_main PUBLIC ${LINK_LIBS} ${TESTU01_LIB} ${TESTU01_MYLIB} ${TESTU01_PROBDIST})
    target_compile_options(testu01_main PUBLIC ${COMPILE_OPTIONS})
endif()

enable_testing()
add_executable(test_rng_loyalty test_rng_loyalty.cc)
target_link_libraries(test_rng_loyalty PUBLIC ${LINK_LIBS})
target_compile_options(test_rng_loyalty PUBLIC ${COMPILE_OPTIONS})
add_test(NAME test_rng_loyalty COMMAND test_rng_loyalty)
