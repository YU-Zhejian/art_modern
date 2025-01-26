cmake_minimum_required(VERSION 3.10)
project(
    htslib
    VERSION 1.21
    LANGUAGES C)

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

find_package(Threads REQUIRED)

include("${LIBCEU_CMAKE_COLLECTION_DIR}/enable_debug.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_c11.cmake")

set(CEU_CM_SILENT ON)
set(CEU_CM_CC_FLAGS "")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/get_compiler_flags.cmake")
set(HTSLIB_CC_FLAGS ${CEU_CM_CC_FLAGS})
unset(CEU_CM_SILENT)

include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_c_helloworld.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_libm.cmake")
# FIXME include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_libpthread.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_libz.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_libbz2.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_liblzma.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/test_libdeflate.cmake")

include_directories("${CMAKE_CURRENT_LIST_DIR}")
include_directories("${CMAKE_CURRENT_LIST_DIR}/htslib")
include_directories("${CMAKE_CURRENT_LIST_DIR}/htscodecs/htscodecs")
include_directories("${CMAKE_BINARY_DIR}/include")
set(LABW_SLIM_HTSLIB_LINK_LIBRARIES "")

if(WIN32
   OR CYGWIN
   OR MSVC
   OR MSYS
   OR MINGW)

else()
    set(HTSLIB_CC_FLAGS -DHAVE_DRAND48 ${HTSLIB_CC_FLAGS})
endif()

if(BUILD_SHARED_LIBS)
    message(NOTICE "Wanting shared htslib.")
    if(NOT ${CEU_CM_HAVE_WORKING_C_HELLOWORLD_RUN_SHARED} EQUAL 0)
        message(FATAL_ERROR "CMake cannot compile & run a hello world program!")
    endif()
    if(NOT ${CEU_CM_HAVE_WORKING_C_C11_RUN_SHARED} EQUAL 0)
        message(FATAL_ERROR "CMake cannot compile & run a C11 program!")
    endif()

    if(${CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_SHARED} EQUAL 0)
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::libm_shared)
    elseif(${CEU_CM_HAVE_WORKING_C_NO_LIBM_RUN_SHARED} EQUAL 0)
        # Do nothing!
    else()
        message(FATAL_ERROR "-lm not found!")
    endif()

    if(NOT ${CEU_CM_HAVE_WORKING_LIBZ_RUN_SHARED} EQUAL 0)
        message(FATAL_ERROR "-lz not found!")
    else()
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::libz_shared)
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBBZ2_RUN_SHARED} EQUAL 0)
        set(HTSLIB_CC_FLAGS -DHAVE_LIBBZ2 ${HTSLIB_CC_FLAGS})

        message(NOTICE "libbz2 detected")
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::libbz2_shared)
    else()
        message(NOTICE "libbz2 NOT detected")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBDEFLATE_RUN_SHARED} EQUAL 0)
        message(NOTICE "libdeflate detected")
        set(HTSLIB_CC_FLAGS -DHAVE_LIBDEFLATE ${HTSLIB_CC_FLAGS})
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::libdeflate_shared)
    else()
        message(NOTICE "libdeflate NOT detected")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBLZMA_RUN_SHARED} EQUAL 0)
        message(NOTICE "liblzma detected")
        set(HTSLIB_CC_FLAGS -DHAVE_LIBLZMA ${HTSLIB_CC_FLAGS})
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::liblzma_shared)
    else()
        message(NOTICE "liblzma NOT detected")
    endif()
else()
    message(NOTICE "Wanting static htslib.")
    if(NOT ${CEU_CM_HAVE_WORKING_C_HELLOWORLD_RUN_STATIC} EQUAL 0)
        message(FATAL_ERROR "CMake cannot compile & run a hello world program!")
    endif()
    if(NOT ${CEU_CM_HAVE_WORKING_C_C11_RUN_STATIC} EQUAL 0)
        message(FATAL_ERROR "CMake cannot compile & run a C11 program!")
    endif()

    if(${CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_STATIC} EQUAL 0)
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::libm_static)
    elseif(${CEU_CM_HAVE_WORKING_C_NO_LIBM_RUN_STATIC} EQUAL 0)
        # Do nothing!
    else()
        message(FATAL_ERROR "-lm not found!")
    endif()

    if(NOT ${CEU_CM_HAVE_WORKING_LIBZ_RUN_STATIC} EQUAL 0)
        message(FATAL_ERROR "-lz not found!")
    else()
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::libz_static)
    endif()
    if(${CEU_CM_HAVE_WORKING_LIBBZ2_RUN_STATIC} EQUAL 0)
        message(NOTICE "libbz2 detected")
        set(HTSLIB_CC_FLAGS -DHAVE_LIBBZ2 ${HTSLIB_CC_FLAGS})
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::libbz2_static)
    else()
        message(NOTICE "libbz2 NOT detected")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBDEFLATE_RUN_STATIC} EQUAL 0)
        message(NOTICE "libdeflate detected")
        set(HTSLIB_CC_FLAGS -DHAVE_LIBDEFLATE ${HTSLIB_CC_FLAGS})
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::libdeflate_static)
    else()
        message(NOTICE "libdeflate NOT detected")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBLZMA_RUN_STATIC} EQUAL 0)
        message(NOTICE "liblzma detected")
        set(HTSLIB_CC_FLAGS -DHAVE_LIBLZMA ${HTSLIB_CC_FLAGS})
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} CEU_CM_EFL::liblzma_static)
    else()
        message(NOTICE "liblzma NOT detected")
    endif()
endif()

configure_file("${CMAKE_CURRENT_LIST_DIR}/htslib_config_vars.h.in" "${CMAKE_BINARY_DIR}/include/htslib_config_vars.h")
configure_file("${CMAKE_CURRENT_LIST_DIR}/htslib_version.h.in" "${CMAKE_BINARY_DIR}/include/htslib_version.h")

file(GLOB HTSLIB_SOURCES ${CMAKE_CURRENT_LIST_DIR}/*.c ${CMAKE_CURRENT_LIST_DIR}/cram/*.c
     ${CMAKE_CURRENT_LIST_DIR}/htscodecs/htscodecs/*.c)
if(BUILD_SHARED_LIBS)
    add_library(labw_slim_htslib SHARED ${HTSLIB_SOURCES})
else()
    add_library(labw_slim_htslib STATIC ${HTSLIB_SOURCES})
    ceu_cm_set_static_target(labw_slim_htslib)
endif()
target_compile_options(labw_slim_htslib BEFORE PRIVATE "${HTSLIB_CC_FLAGS}")

message(STATUS "HTSLib compiler flags: ${HTSLIB_CC_FLAGS}")
message(STATUS "HTSLib link libs: ${LABW_SLIM_HTSLIB_LINK_LIBRARIES}")

target_link_libraries(labw_slim_htslib ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} Threads::Threads)
