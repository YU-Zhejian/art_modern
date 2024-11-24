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
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CYGWIN)
    set(BUILD_SHARED_LIBS OFF)
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

unset(USE_ASIO_PARALLEL)
unset(USE_NOP_PARALLEL)

if(NOT DEFINED USE_ASIO_PARALLEL)
    set(USE_ASIO_PARALLEL ON)
endif()

add_compile_definitions(USE_ASIO_PARALLEL) # USE_NOP_PARALLEL #  TODO

if(NOT DEFINED USE_RANDOM_GENERATOR)
    set(USE_RANDOM_GENERATOR STL)
endif()

set(WITH_MPI OFF) # TODO
set(WITH_PROTOBUF OFF) # TODO

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
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/mpi.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/rng.cmake")

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
