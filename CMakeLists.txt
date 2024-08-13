cmake_minimum_required(VERSION 3.14)
cmake_policy(SET CMP0167 NEW)

project(art_modern LANGUAGES CXX)

# Installation
include(GNUInstallDirs)
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
    set(Boost_USE_STATIC_RUNTIME ON)
endif()
find_package(
    Boost REQUIRED
    COMPONENTS filesystem regex program_options thread log_setup log
    # signals2 is header-only lockfree is header-only
    OPTIONAL_COMPONENTS unit_test_framework timer)
include_directories(${Boost_INCLUDE_DIRS})

# HTSLib
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/labw_slim_htslib")

# Main
configure_file("${CMAKE_CURRENT_LIST_DIR}/test/test_adaptor.h.in" "${CMAKE_BINARY_DIR}/include/test_adaptor.h")

include_directories("${CMAKE_CURRENT_LIST_DIR}/lib")
include_directories("${CMAKE_BINARY_DIR}/include/")
include_directories("${CMAKE_CURRENT_LIST_DIR}/deps/labw_slim_htslib")
file(GLOB ART_MODERN_LIB_SOURCE "lib/*.cc" "lib/*/*.cc")
add_library(art_modern_lib ${ART_MODERN_LIB_SOURCE} lib/art_modern_constants.hh)

set(ART_MODERN_LINK_LIBS
    Boost::filesystem
    Boost::regex
    Boost::program_options
    Boost::thread
    Boost::log_setup
    Boost::log
    hts)
if(Boost_timer_FOUND)
    add_compile_definitions("WITH_BOOST_TIMER")
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::timer)
endif()

file(GLOB ART_SOURCE "art/*.cpp" "art/*.cc")
add_executable(art_modern ${ART_SOURCE})

target_link_libraries(art_modern_lib PUBLIC ${ART_MODERN_LINK_LIBS})

target_link_libraries(art_modern PUBLIC art_modern_lib)

if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_modern_lib)
    ceu_cm_set_static_target(art_modern)
endif()

if(Boost_unit_test_framework_FOUND)
    message(NOTICE "Testing enabled.")
    enable_testing()
    file(GLOB ART_MODERN_TEST_SOURCES "test/*.cc")

    foreach(ART_MODERN_TEST_SOURCE ${ART_MODERN_TEST_SOURCES})
        get_filename_component(ART_MODERN_TEST_EXEC "${ART_MODERN_TEST_SOURCE}" NAME_WLE)
        set(ART_MODERN_TEST_EXEC "art_modern_test_${ART_MODERN_TEST_EXEC}")
        add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
        target_link_libraries("${ART_MODERN_TEST_EXEC}" art_modern_lib)
        add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
    endforeach()
else()
    message(NOTICE "Testing disabled since Boost::unit_test_framework was not found")
endif()
