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
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(CYGWIN)
    set(BUILD_SHARED_LIBS OFF)
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

set(USE_RANDOM_GENERATOR BOOST) # TODO

add_compile_definitions(PARALLEL_DISABLED) # PARALLEL_ASIO TODO

if(NOT DEFINED USE_RANDOM_GENERATOR)
    set(USE_RANDOM_GENERATOR STL)
endif()

unset(USE_STL_RANDOM)
unset(USE_GSL_RANDOM)
unset(USE_BOOST_RANDOM)
unset(USE_ONEMKL_RANDOM)

unset(WITH_MPI)
find_package(MPI REQUIRED) # TODO
include_directories(${MPI_C_HEADER_DIR} ${MPI_CXX_HEADER_DIR})
set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} MPI::MPI_C MPI::MPI_CXX)
set(WITH_MPI OFF) # TODO
set(WITH_PROTOBUF ON) # TODO

if("${USE_RANDOM_GENERATOR}" STREQUAL "STL")
    set(USE_STL_RANDOM ON)
elseif("${USE_RANDOM_GENERATOR}" STREQUAL "GSL")
    set(USE_GSL_RANDOM ON)
    find_package(GSL REQUIRED)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} GSL::gsl)
    include_directories(${GSL_INCLUDE_DIRS})
elseif("${USE_RANDOM_GENERATOR}" STREQUAL "BOOST")
    set(USE_BOOST_RANDOM ON)
elseif("${USE_RANDOM_GENERATOR}" STREQUAL "ONEMKL")
    set(USE_ONEMKL_RANDOM ON)
    find_package(MKL)
    include_directories(${MKL_INCLUDE} ${SYCL_INCLUDE_DIR})
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} MKL::MKL)
else()
    message(
        FATAL_ERROR
            "Variable USE_RANDOM_GENERATOR should be set to one of STL, GSL, BOOST or ONEMKL. Current: ${USE_RANDOM_GENERATOR}"
    )
endif()

# Installation
include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# Macros from LibCEU.
set(LIBCEU_CMAKE_COLLECTION_DIR "${CMAKE_CURRENT_LIST_DIR}/deps/cmake_collections")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/enable_debug.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/libcmake/enhanced_find.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/boost.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/htslib.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/protobuf.cmake")

message(NOTICE "ART_MODERN_LINK_LIBS=${ART_MODERN_LINK_LIBS}")

# Main
configure_file("${CMAKE_CURRENT_LIST_DIR}/test/test_adaptor.h.in" "${CMAKE_BINARY_DIR}/include/test_adaptor.h")
configure_file("${CMAKE_CURRENT_LIST_DIR}/lib/art_modern_config.h.in" "${CMAKE_BINARY_DIR}/include/art_modern_config.h")

include_directories("${CMAKE_CURRENT_LIST_DIR}/lib")
include_directories("${CMAKE_BINARY_DIR}/include/")
include_directories(${env.INCLUDE})

file(GLOB ART_MODERN_LIB_SOURCE "lib/*.cc" "lib/*/*.cc")
add_library(art_modern_lib ${ART_MODERN_LIB_SOURCE} lib/art_modern_constants.hh)

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
