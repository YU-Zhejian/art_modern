cmake_minimum_required(VERSION 3.17)

if(DEFINED USE_CCACHE)
    message(WARNING "-DUSE_CCACHE is deprecated. Ignored.")
endif()

if(DEFINED REPRODUCIBLE_BUILD)
    set(CEU_REPRODUCIBLE_BUILDS
        ON
        CACHE BOOL "Whether to build reproducible builds" FORCE)
else()
    unset(CEU_REPRODUCIBLE_BUILDS)
endif()

project(
    art_modern
    LANGUAGES C CXX
    VERSION 1.1.8
    DESCRIPTION "Modernized ART Simulator of Diverse Next-Generation Sequencing Reads"
    HOMEPAGE_URL "https://github.com/YU-Zhejian/art_modern/")
set(ART_MODERN_VERSION "${PROJECT_VERSION}")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CYGWIN
   OR MINGW
   OR MSYS)
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

set(WITH_MPI OFF) # TODO
set(ART_MODERN_LINK_LIBS "")

# Installation
include(GNUInstallDirs)
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX
        "/usr/local"
        CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# Macros from LibCEU.
set(LIBCEU_CMAKE_COLLECTION_DIR "${CMAKE_CURRENT_LIST_DIR}/deps/cmake_collections")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/enable_debug.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/libcmake/set_static_target.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_libz.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_cxx17.cmake")
unset(WITH_GIT)
include("${LIBCEU_CMAKE_COLLECTION_DIR}/git_adapter.cmake")
if(GIT_FOUND)
    set(WITH_GIT ON)
endif()

# TODO
#if(DEFINED CMAKE_C_FLAGS AND DEFINED CMAKE_CXX_FLAGS)
#    string(REPLACE " " ";" AM_CC_FLAGS "${CMAKE_CXX_FLAGS}")
#else()
unset(CEU_CM_SILENT)
set(CEU_CM_CC_FLAGS "")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/get_compiler_flags.cmake")
set(AM_CC_FLAGS ${CEU_CM_CC_FLAGS})
#endif()

include("${LIBCEU_CMAKE_COLLECTION_DIR}/libcmake/enhanced_find.cmake")

# -----------------------------------------------------------------------------
# Dealing with dependency scripts.
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/boost.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/parallel.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/htslib.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/libfmt.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/absl.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/concurrentqueue.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/qualgen.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/mpi.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/rng.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/malloc.cmake")
# OpenMP temporarily disabled.
# include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/openmp.cmake")

# -----------------------------------------------------------------------------
# Dealing with dependency libs as subdirs.
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/slim_libceu")

# -----------------------------------------------------------------------------
# Dealing with dependency includes.
include_directories(BEFORE "${CMAKE_CURRENT_LIST_DIR}/deps/slim_libceu")
# Threads is now required since lots of
# surveillance and reporting routines are implemented using threads.
find_package(Threads REQUIRED)

if(BUILD_SHARED_LIBS)
    if(DEFINED CEU_CM_HAVE_WORKING_LIBZ_RUN_SHARED AND ${CEU_CM_HAVE_WORKING_LIBZ_RUN_SHARED} EQUAL 0)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} CEU_CM_EFL::libz_shared)
    else()
        message(FATAL_ERROR "-lz not found!")
    endif()
else()
    if(DEFINED CEU_CM_HAVE_WORKING_LIBZ_RUN_STATIC AND ${CEU_CM_HAVE_WORKING_LIBZ_RUN_STATIC} EQUAL 0)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} CEU_CM_EFL::libz_static)
    else()
        message(FATAL_ERROR "-lz not found!")
    endif()
endif()

set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} slim_libceu Threads::Threads)
message(STATUS "art_modern link libs: ${ART_MODERN_LINK_LIBS}")
message(STATUS "art_modern compiler flags: ${AM_CC_FLAGS}")

# -----------------------------------------------------------------------------
# Dealing with config files.

if(NOT EXISTS "${CMAKE_BINARY_DIR}/gen/src/art/builtin_profiles.cc")
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/gen_builtin_profiles.py"
                "${CMAKE_BINARY_DIR}/gen/src/art/"
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" COMMAND_ECHO STDERR
        RESULT_VARIABLE PYTHON_RES)
    if(NOT "${PYTHON_RES}" EQUAL 0)
        message(FATAL_ERROR "Generation of builtin profiles failed!")
    endif()
endif()

if(CEU_CM_SHOULD_ENABLE_TEST)
    configure_file("${CMAKE_CURRENT_LIST_DIR}/src/test/test_adaptor.h.in" "${CMAKE_BINARY_DIR}/include/test_adaptor.h")
endif()
configure_file("${CMAKE_CURRENT_LIST_DIR}/src/libam_support/art_modern_config.h.in"
               "${CMAKE_BINARY_DIR}/include/art_modern_config.h")

# -----------------------------------------------------------------------------
# Dealing with includes.

include_directories(BEFORE "${CMAKE_CURRENT_LIST_DIR}/src/")
include_directories(BEFORE "${CMAKE_BINARY_DIR}/include/")
include_directories(BEFORE "${CMAKE_BINARY_DIR}/gen/src/")

if(DEFINED C_INCLUDE_PATH)
    include_directories(AFTER ${C_INCLUDE_PATH})
endif()

# -----------------------------------------------------------------------------
# Dealing with am_support_lib.

file(GLOB AM_SUPPORT_LIB_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*.cc"
     "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*/*/*.cc" "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*/*.cc")
add_library(am_support_lib ${AM_SUPPORT_LIB_SOURCE})
target_compile_options(am_support_lib BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(am_support_lib PUBLIC ${ART_MODERN_LINK_LIBS})
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(am_support_lib)
endif()
install(TARGETS am_support_lib)

# -----------------------------------------------------------------------------
# Dealing with art_modern executable.

file(GLOB ART_MODERN_LIB_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art/lib/*.cc" "${CMAKE_BINARY_DIR}/gen/src/art/*.cc")
add_library(art_modern_lib ${ART_MODERN_LIB_SOURCE})
target_compile_options(art_modern_lib BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(art_modern_lib PUBLIC am_support_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_modern_lib)
endif()
install(TARGETS art_modern_lib)

file(GLOB ART_MODERN_EXE_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art/exe/*.cc")
add_executable(art_modern ${ART_MODERN_EXE_SOURCE})
target_compile_options(art_modern BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(art_modern PUBLIC art_modern_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_modern)
endif()
install(TARGETS art_modern)

# -----------------------------------------------------------------------------
# Dealing with art_modern tests.

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
        target_link_libraries("${ART_MODERN_TEST_EXEC}" PUBLIC am_support_lib Boost::unit_test_framework)
        add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
        target_compile_options("${ART_MODERN_TEST_EXEC}" BEFORE PRIVATE ${AM_CC_FLAGS})
    endforeach()

    file(GLOB ART_MODERN_TEST_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/art/test/*.cc")

    foreach(ART_MODERN_TEST_SOURCE ${ART_MODERN_TEST_SOURCES})
        get_filename_component(ART_MODERN_TEST_EXEC "${ART_MODERN_TEST_SOURCE}" NAME_WE)
        set(ART_MODERN_TEST_EXEC "am_test_${ART_MODERN_TEST_EXEC}")
        add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
        target_link_libraries("${ART_MODERN_TEST_EXEC}" PUBLIC am_support_lib Boost::unit_test_framework art_modern_lib)
        add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
        target_compile_options("${ART_MODERN_TEST_EXEC}" BEFORE PRIVATE ${AM_CC_FLAGS})
    endforeach()
else()
    message(NOTICE "Testing disabled since Boost::unit_test_framework was not found")
endif()

# -----------------------------------------------------------------------------
# Dealing with art_modern benchmarks.

if(BUILD_ART_MODERN_BENCHMARKS)
    file(GLOB ART_MODERN_BENCHMARK_SOURCES "src/benchmark/*.cc")

    foreach(ART_MODERN_BENCHMARK_SOURCE ${ART_MODERN_BENCHMARK_SOURCES})
        get_filename_component(ART_MODERN_BENCHMARK_EXEC "${ART_MODERN_BENCHMARK_SOURCE}" NAME_WE)
        set(ART_MODERN_BENCHMARK_EXEC "am_bench_${ART_MODERN_BENCHMARK_EXEC}")
        add_executable("${ART_MODERN_BENCHMARK_EXEC}" "${ART_MODERN_BENCHMARK_SOURCE}")
        target_link_libraries("${ART_MODERN_BENCHMARK_EXEC}" PUBLIC am_support_lib)
        target_compile_options("${ART_MODERN_BENCHMARK_EXEC}" BEFORE PRIVATE ${AM_CC_FLAGS})
    endforeach()
endif()

# file(GLOB ART_MODERN_UTIL_SOURCES "src/util/*.cc")

# foreach(ART_MODERN_UTIL_SOURCE ${ART_MODERN_UTIL_SOURCES})
# get_filename_component(ART_MODERN_UTIL_EXEC "${ART_MODERN_UTIL_SOURCE}" NAME_WE)
# set(ART_MODERN_UTIL_EXEC "am_util_${ART_MODERN_UTIL_EXEC}")
# add_executable("${ART_MODERN_UTIL_EXEC}" "${ART_MODERN_UTIL_SOURCE}")
# target_link_libraries("${ART_MODERN_UTIL_EXEC}" PUBLIC am_support_lib)
# endforeach()

# include(CPack)
