cmake_minimum_required(VERSION 3.14)

project(art_modern LANGUAGES CXX)

# Installation
include(GNUInstallDirs)
include(ExternalProject)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# Set C/CXX standard

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# LibCEU
if(DEFINED LIBCEU_SOURCE_PATH AND "${LIBCEU_SOURCE_PATH}" STREQUAL "GET")
    execute_process(
            COMMAND
            git
            clone
            https://github.com/YU-Zhejian/libceu
            "${CMAKE_CURRENT_LIST_DIR}/__libceu"
    )
    set(LIBCEU_SOURCE_PATH "${CMAKE_CURRENT_LIST_DIR}/__libceu")
endif()
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/deps/libceu/" AND NOT DEFINED LIBCEU_SOURCE_PATH)
    set(LIBCEU_SOURCE_PATH "${CMAKE_CURRENT_LIST_DIR}/deps/libceu/")
endif()
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/../libceu/" AND NOT DEFINED LIBCEU_SOURCE_PATH)
    set(LIBCEU_SOURCE_PATH "${CMAKE_CURRENT_LIST_DIR}/../libceu/")
endif()
if(NOT EXISTS "${LIBCEU_SOURCE_PATH}")
    message(FATAL_ERROR "Inclusion of libceu failed.")
endif()

message(STATUS "Inclusion of libceu at ${LIBCEU_SOURCE_PATH}")

ExternalProject_Add(
        "libceu"
        PREFIX "${CMAKE_BINARY_DIR}/__libceu"
        SOURCE_DIR "${LIBCEU_SOURCE_PATH}"
        BINARY_DIR "${CMAKE_BINARY_DIR}/__libceu_build"
        INSTALL_DIR "${CMAKE_BINARY_DIR}/__libceu_install"
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/__libceu_install"
)

include_directories("${CMAKE_BINARY_DIR}/__libceu_install/include")
link_directories("${CMAKE_BINARY_DIR}/__libceu_install/lib")
# add_subdirectory("${LIBCEU_SOURCE_PATH}" __libceu_build)

include("${LIBCEU_SOURCE_PATH}/cmake/libceu.cmake")
add_compile_definitions("HAVE_LIBCEU")

if(CMAKE_VERSION GREATER_EQUAL 3.12)
    add_compile_definitions(CEU_CM_UNDER_CMAKE)
else()
    add_compile_options("-DCEU_CM_UNDER_CMAKE")
endif()

include("${LIBCEU_SOURCE_PATH}/cmake_collections/enable_debug.cmake")

set(BUILD_SHARED_LIBS OFF) # FIXME: Temporary solution under Linux Mint
# Boost
if (NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif ()

if(BUILD_SHARED_LIBS)
    add_definitions(-DBOOST_ALL_DYN_LINK)
    add_definitions(-DBOOST_LOG_DYN_LINK)
    add_definitions(-DBOOST_TEST_DYN_LINK)
else()
    set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(
        Boost
        REQUIRED
        COMPONENTS
        system
        filesystem
        regex
        timer
        program_options
        thread
        log_setup
        log
        unit_test_framework
)
include_directories(${Boost_INCLUDE_DIRS})

# HTSLib
#add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/labw_slim_htslib")
#include("${CMAKE_CURRENT_LIST_DIR}/deps/labw_slim_htslib/cmake/htslib.cmake")

# Main
include_directories("lib")
file(GLOB ART_MODERN_LIB_SOURCE "lib/*.cc")
add_library(art_modern_lib ${ART_MODERN_LIB_SOURCE}
        lib/art_modern_constants.hh)


add_dependencies(art_modern_lib "libceu")
set(
        ART_MODERN_LINK_LIBS
        ${Boost_LIBRARIES}
        art_modern_lib
)
if(BUILD_SHARED_LIBS)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} ${CEU_LIBRARY_SHARED_CXX})
else()
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} ${CEU_LIBRARY_STATIC_CXX})
endif()

file(GLOB ART_SOURCE "art/*.cpp" "art/*.cc")
add_executable(art_modern ${ART_SOURCE})

target_link_libraries(
        art_modern
        PUBLIC
        ${ART_MODERN_LINK_LIBS}
)
enable_testing()


if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_modern_lib)
    ceu_cm_set_static_target(art_modern)
endif()

file(GLOB ART_MODERN_TEST_SOURCES "test/*.cc")

foreach (ART_MODERN_TEST_SOURCE ${ART_MODERN_TEST_SOURCES})
    get_filename_component(ART_MODERN_TEST_EXEC "${ART_MODERN_TEST_SOURCE}" NAME_WLE)
    set(ART_MODERN_TEST_EXEC "art_modern_test_${ART_MODERN_TEST_EXEC}")
    add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
    target_link_libraries("${ART_MODERN_TEST_EXEC}" ${ART_MODERN_LINK_LIBS})
    add_test(COMMAND "${ART_MODERN_TEST_EXEC}" )
endforeach ()
