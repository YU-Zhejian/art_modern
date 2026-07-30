cmake_minimum_required(VERSION 3.17)
set(CEU_VERSION "2026.07.29")
project(
    slim_libceu
    LANGUAGES C CXX
    VERSION ${CEU_VERSION})
include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories("${CMAKE_BINARY_DIR}/include")

# Set C++ standard to C++11 if not already defined
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()
# Do not set C standard.

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)

function(ceu_cm_set_static_target name)
    set_target_properties("${name}" PROPERTIES LINK_SEARCH_START_STATIC 1 LINK_SEARCH_END_STATIC 1 INSTALL_RPATH "")
    if(NOT (BORLAND OR MSVC))
        # NOTE: $<COMPILE_LANGUAGE:...> must not be used here. It is not evaluated per source language in a link
        # context, so it silently dropped -static-libgcc from every C++-linked target. $<LINK_LANGUAGE:...> would be
        # the precise tool, but it requires CMake >= 3.18. Both flags are safe to pass unconditionally instead:
        # -static-libgcc is accepted by the C and the C++ link driver alike, and -static-libstdc++ is ignored by the
        # C driver.
        target_link_options("${name}" PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
endfunction()

if(NOT DEFINED BUILD_SHARED_LIBS)
    if(WIN32
       OR CYGWIN
       OR MINGW
       OR MSVC)
        set(BUILD_SHARED_LIBS OFF)
    else()
        set(BUILD_SHARED_LIBS ON)
    endif()
    message(STATUS "BUILD_SHARED_LIBS is defined to default as: ${BUILD_SHARED_LIBS}")
else()
    message(STATUS "BUILD_SHARED_LIBS is defined externally as: ${BUILD_SHARED_LIBS}")
endif()

if(BUILD_SHARED_LIBS)
    message(STATUS "Building slim_libceu as a shared library")
else()
    message(STATUS "Building slim_libceu as a static library")
endif()

check_include_file("cygwin/version.h" CEU_HAVE_INCLUDE_CYGWIN_VERSION_H)
check_include_file("sys/utsname.h" CEU_HAVE_INCLUDE_SYS_UTSNAME_H)
check_include_file("sys/cdefs.h" CEU_HAVE_INCLUDE_SYS_CDEFS_H)
check_include_file("sys/param.h" CEU_HAVE_INCLUDE_SYS_PARAM_H)
check_include_file("_mingw.h" CEU_HAVE_INCLUDE__MINGW_H)
check_include_file("features.h" CEU_HAVE_INCLUDE_FEATURES_H)
check_include_file("endian.h" CEU_HAVE_INCLUDE_ENDIAN_H)
if(MSVC)
    check_include_file("intrin.h" CEU_HAVE_INCLUDE_INTRIN_H)
else()
    check_include_file("x86intrin.h" CEU_HAVE_INCLUDE_X86INTRIN_H)
endif()

configure_file("${CMAKE_CURRENT_LIST_DIR}/libceu_stddef.h.in" "${CMAKE_BINARY_DIR}/include/libceu_stddef.h")

file(GLOB CEU_LIBRARY_SOURCES ${CMAKE_CURRENT_LIST_DIR}/ceu_check/*.cc ${CMAKE_CURRENT_LIST_DIR}/ceu_check/*.c)
if(BUILD_SHARED_LIBS)
    add_library(slim_libceu SHARED ${CEU_LIBRARY_SOURCES})
else()
    add_library(slim_libceu STATIC ${CEU_LIBRARY_SOURCES})
    ceu_cm_set_static_target(slim_libceu)
endif()
set_target_properties(slim_libceu PROPERTIES VERSION "${PROJECT_VERSION}" SOVERSION 1)
install(TARGETS slim_libceu)

if(DEFINED PROJECT_IS_TOP_LEVEL AND PROJECT_IS_TOP_LEVEL)
    add_executable(ceu_check_exe "${CMAKE_CURRENT_LIST_DIR}/exe/ceu_check_exe.cc")
    target_link_libraries(ceu_check_exe PRIVATE slim_libceu)
    if(NOT BUILD_SHARED_LIBS)
        ceu_cm_set_static_target(ceu_check_exe)
    endif()
else()
    message(STATUS "Skipping ceu_check_exe target because this may not be the top-level project")
endif()
