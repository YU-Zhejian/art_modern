cmake_minimum_required(VERSION 3.14)

project(art_modern LANGUAGES CXX)

# Installation
include(GNUInstallDirs)
include(ExternalProject)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# Set C/CXX standard

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Macros from LibCEU.
macro(ceu_cm_set_static_target name)
    set_target_properties("${name}" PROPERTIES LINK_SEARCH_START_STATIC 1)
    set_target_properties("${name}" PROPERTIES LINK_SEARCH_END_STATIC 1)
    set_target_properties("${name}" PROPERTIES INSTALL_RPATH "")
    if(CMAKE_VERSION GREATER_EQUAL 3.13
       AND NOT BORLAND
       AND NOT MSVC)
        target_link_options(
            "${name}" PRIVATE -static $<$<COMPILE_LANGUAGE:C,CXX>:-static-libgcc>
            $<$<COMPILE_LANGUAGE:CXX>:-static-libstdc++> $<$<COMPILE_LANGUAGE:Fortran>:-static-libgfortran>)
    endif()
endmacro()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

# Boost
if(BUILD_SHARED_LIBS)
    add_definitions(-DBOOST_ALL_DYN_LINK)
    add_definitions(-DBOOST_LOG_DYN_LINK)
    add_definitions(-DBOOST_TEST_DYN_LINK)
else()
    set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(
    Boost REQUIRED
    COMPONENTS system
               filesystem
               regex
               timer
               program_options
               thread
               log_setup
               log
               unit_test_framework)
include_directories(${Boost_INCLUDE_DIRS})

# HTSLib
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/labw_slim_htslib")

# Main
include_directories("lib")
include_directories("${CMAKE_CURRENT_LIST_DIR}/deps/labw_slim_htslib")
file(GLOB ART_MODERN_LIB_SOURCE "lib/*.cc")
add_library(art_modern_lib ${ART_MODERN_LIB_SOURCE} lib/art_modern_constants.hh)

set(ART_MODERN_LINK_LIBS ${Boost_LIBRARIES} hts)

file(GLOB ART_SOURCE "art/*.cpp" "art/*.cc")
add_executable(art_modern ${ART_SOURCE})

target_link_libraries(art_modern_lib PUBLIC ${ART_MODERN_LINK_LIBS})

target_link_libraries(art_modern PUBLIC art_modern_lib)
enable_testing()

if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_modern_lib)
    ceu_cm_set_static_target(art_modern)
endif()

file(GLOB ART_MODERN_TEST_SOURCES "test/*.cc")

foreach(ART_MODERN_TEST_SOURCE ${ART_MODERN_TEST_SOURCES})
    get_filename_component(ART_MODERN_TEST_EXEC "${ART_MODERN_TEST_SOURCE}" NAME_WLE)
    set(ART_MODERN_TEST_EXEC "art_modern_test_${ART_MODERN_TEST_EXEC}")
    add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
    target_link_libraries("${ART_MODERN_TEST_EXEC}" art_modern_lib)
    add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
endforeach()
