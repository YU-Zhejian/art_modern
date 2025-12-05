cmake_minimum_required(VERSION 3.17)

set(DEPRECATED_OPTS USE_CCACHE USE_ABSL USE_BTREE_MAP USE_QUAL_GEN)
foreach(opt ${DEPRECATED_OPTS})
    if(DEFINED ${opt})
        message(WARNING "-D${opt} is deprecated. Ignored.")
        unset(${opt})
    endif()
endforeach()
unset(DEPRECATED_OPTS)

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
    VERSION 1.3.1
    DESCRIPTION "Modernized ART Simulator of Diverse Next-Generation Sequencing Reads"
    HOMEPAGE_URL "https://github.com/YU-Zhejian/art_modern/")
set(ART_MODERN_VERSION "${PROJECT_VERSION}")
set(ART_MODERN_URL "${PROJECT_HOMEPAGE_URL}")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF) # To make C standard conformance stricter.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # To make C++ standard conformance stricter.

if(NOT DEFINED Python3_EXECUTABLE)
    set(Python3_EXECUTABLE python3)
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

if(NOT DEFINED CEU_CM_SHOULD_USE_NATIVE)
    set(CEU_CM_SHOULD_USE_NATIVE OFF)
endif()

if(NOT DEFINED AM_NO_Q_REVERSE)
    set(AM_NO_Q_REVERSE OFF)
endif()

if(NOT DEFINED WITH_MPI)
    set(WITH_MPI OFF)
endif()
set(ART_MODERN_LINK_LIBS "")

# Installation
include(GNUInstallDirs)
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
    message(WARNING "CMAKE_INSTALL_PREFIX not defined, setting to /usr/local")
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
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_cxx_stdthread.cmake")
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

# -----------------------------------------------------------------------------
# Dealing with dependency scripts.
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/boost.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/parallel.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/htslib.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/libfmt.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/concurrentqueue.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/mpi.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/rng.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/src/cmake.d/malloc.cmake")

# -----------------------------------------------------------------------------
# Dealing with dependency libs as subdirs.
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/slim_libceu")

# -----------------------------------------------------------------------------
# Dealing with dependency includes.
include_directories(BEFORE "${CMAKE_CURRENT_LIST_DIR}/deps/slim_libceu")
include_directories(BEFORE "${CMAKE_CURRENT_LIST_DIR}/deps/cgranges/cpp")
include_directories(BEFORE "${CMAKE_CURRENT_LIST_DIR}/deps/EWAHBoolArray/include")

if(BUILD_SHARED_LIBS)
    if(DEFINED CEU_CM_HAVE_WORKING_LIBZ_RUN_SHARED AND ${CEU_CM_HAVE_WORKING_LIBZ_RUN_SHARED} EQUAL 0)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} CEU_CM_EFL::libz_shared)
    else()
        message(FATAL_ERROR "-lz not found!")
    endif()
    if(DEFINED CEU_CM_HAVE_WORKING_CXX_STDTHREAD_RUN_SHARED AND ${CEU_CM_HAVE_WORKING_CXX_STDTHREAD_RUN_SHARED} EQUAL 0)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Threads::Threads)
    else()
        message(FATAL_ERROR "CXX stdthread not found!")
    endif()
else()
    if(DEFINED CEU_CM_HAVE_WORKING_LIBZ_RUN_STATIC AND ${CEU_CM_HAVE_WORKING_LIBZ_RUN_STATIC} EQUAL 0)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} CEU_CM_EFL::libz_static)
    else()
        message(FATAL_ERROR "-lz not found!")
    endif()
    if(DEFINED CEU_CM_HAVE_WORKING_CXX_STDTHREAD_RUN_STATIC AND ${CEU_CM_HAVE_WORKING_CXX_STDTHREAD_RUN_STATIC} EQUAL 0)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Threads::Threads)
    else()
        message(FATAL_ERROR "CXX stdthread not found!")
    endif()
endif()

set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} slim_libceu)
message(STATUS "art_modern link libs: ${ART_MODERN_LINK_LIBS}")
message(STATUS "art_modern compiler flags: ${AM_CC_FLAGS}")

# -----------------------------------------------------------------------------
# Dealing with config files.

if(NOT EXISTS "${CMAKE_BINARY_DIR}/gen/src/art/builtin_profiles.c")
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
# Dealing with ${AM_SUPPORT_LIB}.
if(WITH_MPI)
    set(AM_MPI_POSTFIX "-mpi")
else()
    set(AM_MPI_POSTFIX "")
endif()

set(AM_SUPPORT_LIB "am_support_lib${AM_MPI_POSTFIX}")
file(
    GLOB
    AM_SUPPORT_LIB_SOURCE
    "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*.cc"
    "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*/*/*.cc"
    "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*/*.cc"
    "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*.c"
    "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*/*/*.c"
    "${CMAKE_CURRENT_LIST_DIR}/src/libam_support/*/*.c")
add_library(${AM_SUPPORT_LIB} ${AM_SUPPORT_LIB_SOURCE})
target_compile_options(${AM_SUPPORT_LIB} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${AM_SUPPORT_LIB} PUBLIC ${ART_MODERN_LINK_LIBS})
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${AM_SUPPORT_LIB})
endif()
install(TARGETS ${AM_SUPPORT_LIB})

# -----------------------------------------------------------------------------
# Dealing with art_modern executable.
set(ART_MODERN_LIB "art_modern_lib${AM_MPI_POSTFIX}")
file(GLOB ART_MODERN_LIB_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art/lib/*.cc" "${CMAKE_BINARY_DIR}/gen/src/art/*.c")
add_library(${ART_MODERN_LIB} ${ART_MODERN_LIB_SOURCE})
target_compile_options(${ART_MODERN_LIB} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${ART_MODERN_LIB} PUBLIC ${AM_SUPPORT_LIB})
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${ART_MODERN_LIB})
endif()
install(TARGETS ${ART_MODERN_LIB})

set(ART_MODERN_EXE "art_modern${AM_MPI_POSTFIX}")
file(GLOB ART_MODERN_EXE_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art/exe/*.cc")
add_executable(${ART_MODERN_EXE} ${ART_MODERN_EXE_SOURCE})
target_compile_options(${ART_MODERN_EXE} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${ART_MODERN_EXE} PUBLIC ${ART_MODERN_LIB})
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${ART_MODERN_EXE})
endif()
install(TARGETS ${ART_MODERN_EXE})

set(ART_PROFILE_BUILDER_LIB "art_profile_builder_lib${AM_MPI_POSTFIX}")
file(GLOB ART_PROFILE_BUILDER_LIB_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art_profile_builder/lib/*.cc")
add_library(${ART_PROFILE_BUILDER_LIB} ${ART_PROFILE_BUILDER_LIB_SOURCE})
target_compile_options(${ART_PROFILE_BUILDER_LIB} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${ART_PROFILE_BUILDER_LIB} PUBLIC ${AM_SUPPORT_LIB})
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${ART_PROFILE_BUILDER_LIB})
endif()
install(TARGETS ${ART_PROFILE_BUILDER_LIB})

set(ART_PROFILER_BUILDER_EXE "art_profile_builder${AM_MPI_POSTFIX}")
file(GLOB ART_PROFILE_BUILDER_EXE_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art_profile_builder/exe/*.cc")
add_executable(${ART_PROFILE_BUILDER_EXE} ${ART_PROFILE_BUILDER_EXE_SOURCE})
target_compile_options(${ART_PROFILE_BUILDER_EXE} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${ART_PROFILE_BUILDER_EXE} PUBLIC art_profile_builder_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${ART_PROFILE_BUILDER_EXE})
endif()
install(TARGETS ${ART_PROFILE_BUILDER_EXE})

file(GLOB ART_TSAM2GSAM_LIB_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art_tsam2gsam/lib/*.cc")
add_library(art_tsam2gsam_lib ${ART_TSAM2GSAM_LIB_SOURCE})
target_compile_options(art_tsam2gsam_lib BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(art_tsam2gsam_lib PUBLIC am_support_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_tsam2gsam_lib)
endif()
install(TARGETS art_tsam2gsam_lib)

if(WITH_MPI)
    set(ART_TSAM2GSAM_EXE art_tsam2gsam-mpi)
else()
    set(ART_TSAM2GSAM_EXE art_tsam2gsam)
endif()

file(GLOB ART_TSAM2GSAM_EXE_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art_tsam2gsam/exe/*.cc")
add_executable(${ART_TSAM2GSAM_EXE} ${ART_TSAM2GSAM_EXE_SOURCE})
target_compile_options(${ART_TSAM2GSAM_EXE} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${ART_TSAM2GSAM_EXE} PUBLIC art_tsam2gsam_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${ART_TSAM2GSAM_EXE})
endif()
install(TARGETS ${ART_TSAM2GSAM_EXE})

file(GLOB ART_SAMCMPS_LIB_SOURCE "${CMAKE_CURRENT_LIST_DIR}/src/art_samcmps/lib/*.cc")
add_library(art_samcmps_lib ${ART_SAMCMPS_LIB_SOURCE})
target_compile_options(art_samcmps_lib BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(art_samcmps_lib PUBLIC am_support_lib art_tsam2gsam_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_samcmps_lib)
endif()
install(TARGETS art_samcmps_lib)

if(WITH_MPI)
    set(ART_SAMCMP_EXE art_samcmp-mpi)
else()
    set(ART_SAMCMP_EXE art_samcmp)
endif()

add_executable(${ART_SAMCMP_EXE} "${CMAKE_CURRENT_LIST_DIR}/src/art_samcmps/exe/samcmp.cc")
target_compile_options(${ART_SAMCMP_EXE} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${ART_SAMCMP_EXE} PUBLIC art_samcmps_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${ART_SAMCMP_EXE})
endif()
install(TARGETS ${ART_SAMCMP_EXE})

if(WITH_MPI)
    set(ART_SAMCMP3_EXE art_samcmp3-mpi)
else()
    set(ART_SAMCMP3_EXE art_samcmp3)
endif()

add_executable(${ART_SAMCMP3_EXE} "${CMAKE_CURRENT_LIST_DIR}/src/art_samcmps/exe/samcmp3.cc")
target_compile_options(${ART_SAMCMP3_EXE} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${ART_SAMCMP3_EXE} PUBLIC art_samcmps_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${ART_SAMCMP3_EXE})
endif()
install(TARGETS ${ART_SAMCMP3_EXE})

if(WITH_MPI)
    set(ART_SAMGTFCMP_EXE art_samgtfcmp-mpi)
else()
    set(ART_SAMGTFCMP_EXE art_samgtfcmp)
endif()

add_executable(${ART_SAMGTFCMP_EXE} "${CMAKE_CURRENT_LIST_DIR}/src/art_samcmps/exe/samgtfcmp.cc")
target_compile_options(${ART_SAMGTFCMP_EXE} BEFORE PRIVATE ${AM_CC_FLAGS})
target_link_libraries(${ART_SAMGTFCMP_EXE} PUBLIC art_samcmps_lib)
if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(${ART_SAMGTFCMP_EXE})
endif()
install(TARGETS ${ART_SAMGTFCMP_EXE})

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
        set(ART_MODERN_TEST_EXEC "am_test_${ART_MODERN_TEST_EXEC}${AM_MPI_POSTFIX}")
        add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
        target_link_libraries("${ART_MODERN_TEST_EXEC}" PUBLIC ${AM_SUPPORT_LIB} Boost::unit_test_framework)
        target_compile_options("${ART_MODERN_TEST_EXEC}" BEFORE PRIVATE ${AM_CC_FLAGS})
        if(NOT BUILD_SHARED_LIBS)
            ceu_cm_set_static_target("${ART_MODERN_TEST_EXEC}")
        endif()
        add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
    endforeach()

    file(GLOB ART_MODERN_TEST_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/art/test/*.cc")
    foreach(ART_MODERN_TEST_SOURCE ${ART_MODERN_TEST_SOURCES})
        get_filename_component(ART_MODERN_TEST_EXEC "${ART_MODERN_TEST_SOURCE}" NAME_WE)
        set(ART_MODERN_TEST_EXEC "am_test_${ART_MODERN_TEST_EXEC}${AM_MPI_POSTFIX}")
        add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
        target_link_libraries("${ART_MODERN_TEST_EXEC}" PUBLIC ${AM_SUPPORT_LIB} Boost::unit_test_framework
                                                               ${ART_MODERN_LIB})
        add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
        if(NOT BUILD_SHARED_LIBS)
            ceu_cm_set_static_target("${ART_MODERN_TEST_EXEC}")
        endif()
        target_compile_options("${ART_MODERN_TEST_EXEC}" BEFORE PRIVATE ${AM_CC_FLAGS})
    endforeach()
else()
    message(NOTICE "Testing disabled since Boost::unit_test_framework was not found")
endif()

# -----------------------------------------------------------------------------
# Dealing with art_modern benchmarks.

if(DEFINED BUILD_ART_MODERN_BENCHMARKS AND BUILD_ART_MODERN_BENCHMARKS)
    configure_file("${CMAKE_CURRENT_LIST_DIR}/src/benchmark/benchmark_adaptor.h.in"
                   "${CMAKE_BINARY_DIR}/include/benchmark_adaptor.h")
    file(GLOB ART_MODERN_BENCHMARK_SOURCES "src/benchmark/*.cc")
    foreach(ART_MODERN_BENCHMARK_SOURCE ${ART_MODERN_BENCHMARK_SOURCES})
        get_filename_component(ART_MODERN_BENCHMARK_EXEC "${ART_MODERN_BENCHMARK_SOURCE}" NAME_WE)
        set(ART_MODERN_BENCHMARK_EXEC "am_bench_${ART_MODERN_BENCHMARK_EXEC}${AM_MPI_POSTFIX}")
        add_executable("${ART_MODERN_BENCHMARK_EXEC}" "${ART_MODERN_BENCHMARK_SOURCE}")
        target_link_libraries("${ART_MODERN_BENCHMARK_EXEC}" PUBLIC ${AM_SUPPORT_LIB})
        target_compile_options("${ART_MODERN_BENCHMARK_EXEC}" BEFORE PRIVATE ${AM_CC_FLAGS})
        if(NOT BUILD_SHARED_LIBS)
            ceu_cm_set_static_target("${ART_MODERN_TEST_EXEC}")
        endif()
    endforeach()
endif()
