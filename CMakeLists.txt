cmake_minimum_required(VERSION 3.14)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
    cmake_policy(SET CMP0167 NEW)
endif()

project(
    art_modern
    LANGUAGES C CXX
    VERSION 1.0.0)
set(ART_MODERN_VERSION "${PROJECT_VERSION}")

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

# Installation
include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# Macros from LibCEU.
set(LIBCEU_CMAKE_COLLECTION_DIR "${CMAKE_CURRENT_LIST_DIR}/deps/cmake_collections")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/enable_debug.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/libcmake/enhanced_find.cmake")

if(Boost_stacktrace_backtrace_FOUND)
    add_definitions(-DBOOST_STACKTRACE_USE_BACKTRACE)
endif()

# Boost
if(BUILD_SHARED_LIBS)
    add_definitions(-DBOOST_ALL_DYN_LINK)
    add_definitions(-DBOOST_LOG_DYN_LINK)
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DBOOST_STACKTRACE_DYN_LINK)
    add_definitions(-DBOOST_STACKTRACE_LINK)
else()
    set(Boost_USE_STATIC_LIBS ON)
    # set(Boost_USE_STATIC_RUNTIME ON)
endif()
find_package(
    Boost REQUIRED
    COMPONENTS filesystem
               regex
               program_options
               thread
               log_setup
               log
    # signals2 is header-only lockfree is header-only
    OPTIONAL_COMPONENTS unit_test_framework timer stacktrace_basic stacktrace_backtrace stacktrace_windbg)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${env.INCLUDE})

# HTSLib
if(DEFINED USE_HTSLIB)

    if(BUILD_SHARED_LIBS)
        ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CMAKE_CM_LIBHTS_LIBRARY LINKER_FLAG hts PKGCONFIG_NAME htslib)
    else()
        ceu_cm_enhanced_find_library(
            OUTPUT_VARIABLE
            CMAKE_CM_LIBHTS_LIBRARY
            STATIC
            LINKER_FLAG
            hts
            PKGCONFIG_NAME
            htslib)
    endif()

    if(NOT CMAKE_CM_LIBHTS_LIBRARY)
        message(FATAL_ERROR "htslib not found!")
    endif()

    set(USED_HTSLIB_NAME "${CMAKE_CM_LIBHTS_LIBRARY}")
else()
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/labw_slim_htslib")
    include_directories("${CMAKE_CURRENT_LIST_DIR}/deps/labw_slim_htslib")
    set(USED_HTSLIB_NAME labw_slim_htslib)
endif()
add_compile_definitions(USE_STD_RANDOM)
add_compile_definitions(PARALLEL_ASIO)

# find_package(MKL) # TODO find_package(IntelSYCL) # TODO include_directories(${MKL_INCLUDE} ${SYCL_INCLUDE_DIR}) # TODO

# Main
configure_file("${CMAKE_CURRENT_LIST_DIR}/test/test_adaptor.h.in" "${CMAKE_BINARY_DIR}/include/test_adaptor.h")
configure_file("${CMAKE_CURRENT_LIST_DIR}/lib/art_modern_config.h.in" "${CMAKE_BINARY_DIR}/include/art_modern_config.h")

include_directories("${CMAKE_CURRENT_LIST_DIR}/lib")
include_directories("${CMAKE_BINARY_DIR}/include/")
file(GLOB ART_MODERN_LIB_SOURCE "lib/*.cc" "lib/*/*.cc")
add_library(art_modern_lib ${ART_MODERN_LIB_SOURCE} lib/art_modern_constants.hh)

set(ART_MODERN_LINK_LIBS
    Boost::filesystem
    Boost::regex
    Boost::program_options
    Boost::thread
    Boost::log_setup
    Boost::log
    "${USED_HTSLIB_NAME}")
if(Boost_timer_FOUND)
    add_compile_definitions("WITH_BOOST_TIMER")
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::timer)
endif()
if(Boost_stacktrace_backtrace_FOUND)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::stacktrace_backtrace)
elseif(Boost_stacktrace_windbg_FOUND)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::stacktrace_windbg)
elseif (Boost_stacktrace_basic_FOUND)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} Boost::stacktrace_basic)
endif()

file(GLOB ART_SOURCE "art/*.cpp" "art/*.cc")
add_executable(art_modern ${ART_SOURCE})

target_link_libraries(art_modern_lib PUBLIC ${ART_MODERN_LINK_LIBS})

target_link_libraries(art_modern PUBLIC art_modern_lib)

if(NOT BUILD_SHARED_LIBS)
    ceu_cm_set_static_target(art_modern_lib)
    ceu_cm_set_static_target(art_modern)
endif()
if(NOT CEU_CM_SHOULD_ENABLE_TEST)
    message(NOTICE "Testing disabled since CEU_CM_SHOULD_ENABLE_TEST was set to FASE")
elseif(Boost_unit_test_framework_FOUND)
    message(NOTICE "Testing enabled.")
    enable_testing()
    file(GLOB ART_MODERN_TEST_SOURCES "test/*.cc")

    foreach(ART_MODERN_TEST_SOURCE ${ART_MODERN_TEST_SOURCES})
        get_filename_component(ART_MODERN_TEST_EXEC "${ART_MODERN_TEST_SOURCE}" NAME_WLE)
        set(ART_MODERN_TEST_EXEC "art_modern_test_${ART_MODERN_TEST_EXEC}")
        add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
        target_link_libraries("${ART_MODERN_TEST_EXEC}" art_modern_lib Boost::unit_test_framework)
        add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
    endforeach()
else()
    message(NOTICE "Testing disabled since Boost::unit_test_framework was not found")
endif()
