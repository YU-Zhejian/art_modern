cmake_minimum_required(VERSION 3.10)

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

if(NOT DEFINED USE_CCACHE)
    set(USE_CCACHE OFF)
endif()
if(USE_CCACHE)
    include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/ccache.cmake")
endif()

project(
    art_modern
    LANGUAGES C CXX
    VERSION 1.1.1)
set(ART_MODERN_VERSION "${PROJECT_VERSION}")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CYGWIN)
    set(BUILD_SHARED_LIBS OFF)
endif()

if(NOT DEFINED Python3_EXECUTABLE)
    set(Python3_EXECUTABLE python3)
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

if(NOT DEFINED CEU_CM_SHOULD_USE_NATIVE)
    set(CEU_CM_SHOULD_USE_NATIVE OFF)
endif()

if(NOT DEFINED USE_THREAD_PARALLEL)
    set(USE_THREAD_PARALLEL ASIO)
endif()

if(NOT DEFINED USE_BTREE_MAP)
    set(USE_BTREE_MAP ON)
endif()

if(NOT DEFINED USE_RANDOM_GENERATOR)
    set(USE_RANDOM_GENERATOR STL)
endif()

if(NOT DEFINED USE_MALLOC)
    set(USE_MALLOC AUTO)
endif()

set(WITH_MPI OFF) # TODO

# Installation
include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# Macros from LibCEU.
set(LIBCEU_CMAKE_COLLECTION_DIR "${CMAKE_CURRENT_LIST_DIR}/deps/cmake_collections")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/enable_debug.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/libcmake/set_static_target.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_cxx17.cmake")
unset(WITH_GIT)
include("${LIBCEU_CMAKE_COLLECTION_DIR}/git_adapter.cmake")
if(GIT_FOUND)
    set(WITH_GIT ON)
endif()

unset(CEU_CM_SILENT)
set(CEU_CM_CC_FLAGS "")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/get_compiler_flags.cmake")
set(AM_CC_FLAGS ${CEU_CM_CC_FLAGS})

include("${LIBCEU_CMAKE_COLLECTION_DIR}/libcmake/enhanced_find.cmake")

set(ART_MODERN_LINK_LIBS "")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/boost.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/parallel.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/htslib.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/mpi.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/rng.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/malloc.cmake")
# OpenMP temporarily disabled.
# include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/openmp.cmake")

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/slim_libceu")
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/slim_fmt")
include_directories(BEFORE SYSTEM "${CMAKE_CURRENT_LIST_DIR}/deps/slim_libceu")
include_directories(BEFORE SYSTEM "${CMAKE_CURRENT_LIST_DIR}/deps/concurrentqueue")
# include_directories("${CMAKE_CURRENT_LIST_DIR}/deps/simde")
include_directories(BEFORE SYSTEM "${CMAKE_CURRENT_LIST_DIR}/deps/base64/")
include_directories(BEFORE SYSTEM "${CMAKE_CURRENT_LIST_DIR}/deps/cpp-btree/")
include_directories(BEFORE SYSTEM "${CMAKE_CURRENT_LIST_DIR}/deps/slim_fmt/include/")
include_directories(BEFORE SYSTEM "${CMAKE_CURRENT_LIST_DIR}/deps/thread-pool/include/")
include_directories(BEFORE SYSTEM "${CMAKE_CURRENT_LIST_DIR}/deps/slim_abseil/")

# Threads is now required since lots of
# surveillance and reporting routines are implemented using threads.
find_package(Threads REQUIRED)

set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} slim_libceu slim_libfmt Threads::Threads)
message(STATUS "art_modern link libs: ${ART_MODERN_LINK_LIBS}")
message(STATUS "art_modern compiler flags: ${AM_CC_FLAGS}")

if(NOT EXISTS "${CMAKE_BINARY_DIR}/gen/src/art/builtin_profiles.cc")
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/sh.d/gen_builtin_profiles.py"
                "${CMAKE_BINARY_DIR}/gen/src/art/"
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" COMMAND_ECHO STDERR
        RESULT_VARIABLE PYTHON_RES)
    if(NOT "${PYTHON_RES}" EQUAL 0)
        message(FATAL_ERROR "Generation of builtin profiles failed!")
    endif()
endif()
# Main
if(CEU_CM_SHOULD_ENABLE_TEST)
    configure_file("${CMAKE_CURRENT_LIST_DIR}/src/test/test_adaptor.h.in" "${CMAKE_BINARY_DIR}/include/test_adaptor.h")
endif()
configure_file("${CMAKE_CURRENT_LIST_DIR}/src/lib/art_modern_config.h.in"
               "${CMAKE_BINARY_DIR}/include/art_modern_config.h")

include_directories(BEFORE "${CMAKE_CURRENT_LIST_DIR}/src/")
include_directories(BEFORE "${CMAKE_BINARY_DIR}/include/")
include_directories(BEFORE "${CMAKE_BINARY_DIR}/gen/src/")

if(DEFINED env.INCLUDE)
    include_directories(AFTER SYSTEM ${env.INCLUDE})
endif()
if(DEFINED CMAKE_INCLUDE_PATH)
    include_directories(AFTER SYSTEM ${CMAKE_INCLUDE_PATH})
endif()

file(GLOB ART_MODERN_LIB_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/libam/*.cc"
     "${CMAKE_CURRENT_LIST_DIR}/src/libam/*/*/*.cc" "${CMAKE_CURRENT_LIST_DIR}/src/libam/*/*.cc")
add_library(art_modern_lib ${ART_MODERN_LIB_SOURCE})
target_compile_options(art_modern_lib BEFORE PRIVATE ${AM_CC_FLAGS})

file(GLOB ART_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art/*.cc" "${CMAKE_BINARY_DIR}/gen/src/art/*.cc")
add_executable(art_modern ${ART_SOURCE})
target_compile_options(art_modern BEFORE PRIVATE ${AM_CC_FLAGS})

target_link_libraries(art_modern_lib PRIVATE ${ART_MODERN_LINK_LIBS})
target_link_libraries(art_modern PRIVATE art_modern_lib)

if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_modern_lib)
    ceu_cm_set_static_target(art_modern)
endif()

if(NOT CEU_CM_SHOULD_ENABLE_TEST)
    message(NOTICE "Testing disabled since CEU_CM_SHOULD_ENABLE_TEST was set to FASE")
elseif(Boost_unit_test_framework_FOUND)
    message(NOTICE "Testing enabled.")
    enable_testing()
    file(GLOB ART_MODERN_TEST_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/test/*.cc")

    foreach(ART_MODERN_TEST_SOURCE ${ART_MODERN_TEST_SOURCES})
        get_filename_component(ART_MODERN_TEST_EXEC "${ART_MODERN_TEST_SOURCE}" NAME_WE)
        set(ART_MODERN_TEST_EXEC "am_test_${ART_MODERN_TEST_EXEC}")
        add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
        target_link_libraries("${ART_MODERN_TEST_EXEC}" art_modern_lib Boost::unit_test_framework)
        add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
        target_compile_options("${ART_MODERN_TEST_EXEC}" BEFORE PRIVATE ${AM_CC_FLAGS})
    endforeach()
else()
    message(NOTICE "Testing disabled since Boost::unit_test_framework was not found")
endif()

if(CEU_CM_SHOULD_ENABLE_TEST)
    file(GLOB ART_MODERN_BENCHMARK_SOURCES "src/benchmark/*.cc")

    foreach(ART_MODERN_BENCHMARK_SOURCE ${ART_MODERN_BENCHMARK_SOURCES})
        get_filename_component(ART_MODERN_BENCHMARK_EXEC "${ART_MODERN_BENCHMARK_SOURCE}" NAME_WE)
        set(ART_MODERN_BENCHMARK_EXEC "am_bench_${ART_MODERN_BENCHMARK_EXEC}")
        add_executable("${ART_MODERN_BENCHMARK_EXEC}" "${ART_MODERN_BENCHMARK_SOURCE}")
        target_link_libraries("${ART_MODERN_BENCHMARK_EXEC}" art_modern_lib)
        target_compile_options("${ART_MODERN_BENCHMARK_EXEC}" BEFORE PRIVATE ${AM_CC_FLAGS})
    endforeach()
endif()

# file(GLOB ART_MODERN_UTIL_SOURCES "src/util/*.cc")

# foreach(ART_MODERN_UTIL_SOURCE ${ART_MODERN_UTIL_SOURCES})
# get_filename_component(ART_MODERN_UTIL_EXEC "${ART_MODERN_UTIL_SOURCE}" NAME_WE)
# set(ART_MODERN_UTIL_EXEC "am_util_${ART_MODERN_UTIL_EXEC}")
# add_executable("${ART_MODERN_UTIL_EXEC}" "${ART_MODERN_UTIL_SOURCE}")
# target_link_libraries("${ART_MODERN_UTIL_EXEC}" art_modern_lib)
# endforeach()

# include(CPack)
