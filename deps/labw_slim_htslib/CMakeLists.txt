cmake_minimum_required(VERSION 3.17)
project(
    htslib
    VERSION 1.21
    LANGUAGES C)

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

find_package(Threads REQUIRED)

include("${LIBCEU_CMAKE_COLLECTION_DIR}/enable_debug.cmake")
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

if(WIN32 OR CYGWIN OR MSVC OR MSYS OR MINGW)
else()
    add_compile_definitions(HAVE_DRAND48)
    # It seems that MinGW also requries -lregex.
endif()

if(BUILD_SHARED_LIBS)
    message(NOTICE "Wanting shared htslib.")
    if(NOT ${CEU_CM_HAVE_WORKING_C_HELLOWORLD_RUN_SHARED} EQUAL 0)
        message(FATAL_ERROR "CMake cannot compile & run a hello world program!")
    endif()

    if(NOT ${CEU_CM_HAVE_WORKING_LIBZ_RUN_SHARED} EQUAL 0)
        message(FATAL_ERROR "-lz not found!")
    else()
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} "${CEU_CM_LIBZ_LIBRARY_SHARED}")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBBZ2_RUN_SHARED} EQUAL 0)
        message(NOTICE "libbz2 detected")
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} "${CEU_CM_LIBBZ2_LIBRARY_SHARED}")
        add_compile_definitions(HAVE_LIBBZ2)
    else()
        message(NOTICE "libbz2 NOT detected")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBDEFLATE_RUN_SHARED} EQUAL 0)
        message(NOTICE "libdeflate detected")
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} "${CEU_CM_LIBDEFLATE_LIBRARY_SHARED}")
        add_compile_definitions(HAVE_LIBDEFLATE)
    else()
        message(NOTICE "libdeflate NOT detected")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBLZMA_RUN_SHARED} EQUAL 0)
        message(NOTICE "liblzma detected")
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} "${CEU_CM_LIBLZMA_LIBRARY_SHARED}")
        add_compile_definitions(HAVE_LIBLZMA)
    else()
        message(NOTICE "liblzma NOT detected")
    endif()
else()
    message(NOTICE "Wanting static htslib.")
    if(NOT ${CEU_CM_HAVE_WORKING_C_HELLOWORLD_RUN_STATIC} EQUAL 0)
        message(FATAL_ERROR "CMake cannot compile & run a hello world program!")
    endif()

    if(NOT ${CEU_CM_HAVE_WORKING_LIBZ_RUN_STATIC} EQUAL 0)
        message(FATAL_ERROR "-lz not found!")
    else()
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} "${CEU_CM_LIBZ_LIBRARY_STATIC}")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBBZ2_RUN_STATIC} EQUAL 0)
        message(NOTICE "libbz2 detected")
        add_compile_definitions(HAVE_LIBBZ2)
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} "${CEU_CM_LIBBZ2_LIBRARY_STATIC}")
    else()
        message(NOTICE "libbz2 NOT detected")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBDEFLATE_RUN_STATIC} EQUAL 0)
        message(NOTICE "libdeflate detected")
        add_compile_definitions(HAVE_LIBDEFLATE)
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} "${CEU_CM_LIBDEFLATE_LIBRARY_STATIC}")
    else()
        message(NOTICE "libdeflate NOT detected")
    endif()

    if(${CEU_CM_HAVE_WORKING_LIBLZMA_RUN_STATIC} EQUAL 0)
        message(NOTICE "liblzma detected")
        add_compile_definitions(HAVE_LIBLZMA)
        set(LABW_SLIM_HTSLIB_LINK_LIBRARIES ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} "${CEU_CM_LIBLZMA_LIBRARY_STATIC}")
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

target_link_libraries(
        labw_slim_htslib
        ${LABW_SLIM_HTSLIB_LINK_LIBRARIES} Threads::Threads)
