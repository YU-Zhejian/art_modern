cmake_minimum_required(VERSION 3.17)
project(
    slim_libceu
    LANGUAGES C CXX
    VERSION 2025.12.03)
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

check_include_file("cygwin/version.h" CEU_HAVE_INCLUDE_CYGWIN_VERSION_H)
check_include_file("sys/utsname.h" CEU_HAVE_INCLUDE_SYS_UTSNAME_H)
check_include_file("sys/cdefs.h" CEU_HAVE_INCLUDE_SYS_CDEFS_H)
check_include_file("_mingw.h" CEU_HAVE_INCLUDE__MINGW_H)
check_include_file("features.h" CEU_HAVE_INCLUDE_FEATURES_H)

configure_file("${CMAKE_CURRENT_LIST_DIR}/libceu_stddef.h.in" "${CMAKE_BINARY_DIR}/include/libceu_stddef.h")

file(GLOB CEU_LIBRARY_SOURCES ${CMAKE_CURRENT_LIST_DIR}/*/*.cc)
if(BUILD_SHARED_LIBS)
    add_library(slim_libceu SHARED ${CEU_LIBRARY_SOURCES})
else()
    add_library(slim_libceu STATIC ${CEU_LIBRARY_SOURCES})
endif()
set_target_properties(slim_libceu PROPERTIES VERSION "${PROJECT_VERSION}" SOVERSION 1)
install(TARGETS slim_libceu)

file(GLOB CEU_LIBRARY_SOURCES ${CMAKE_CURRENT_LIST_DIR}/*/*.c)
if(BUILD_SHARED_LIBS)
    add_library(slim_libceu_c SHARED ${CEU_LIBRARY_SOURCES})
else()
    add_library(slim_libceu_c STATIC ${CEU_LIBRARY_SOURCES})
endif()
set_target_properties(slim_libceu_c PROPERTIES VERSION "${PROJECT_VERSION}" SOVERSION 1)
install(TARGETS slim_libceu_c)

if(DEFINED PROJECT_IS_TOP_LEVEL AND PROJECT_IS_TOP_LEVEL)
    add_executable(ceu_check_exe exe/ceu_check_exe.cc)
    target_link_libraries(ceu_check_exe PUBLIC slim_libceu)
    add_executable(ceu_check_cc_c exe/ceu_check_cc.c)
    target_link_libraries(ceu_check_cc_c PUBLIC slim_libceu_c)
else()
    message(STATUS "Skipping ceu_check_exe target because this may not be the top-level project")
endif()
