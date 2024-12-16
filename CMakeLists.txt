cmake_minimum_required(VERSION 3.17)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
    cmake_policy(SET CMP0167 NEW)
endif()

project(
    art_modern
    LANGUAGES C CXX
    VERSION 1.0.0)
set(ART_MODERN_VERSION "${PROJECT_VERSION}")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CYGWIN)
    set(BUILD_SHARED_LIBS OFF)
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

if(NOT DEFINED USE_THREAD_PARALLEL)
    set(USE_THREAD_PARALLEL ASIO)
endif()

if(NOT DEFINED USE_RANDOM_GENERATOR)
    set(USE_RANDOM_GENERATOR STL)
endif()

set(WITH_MPI OFF) # TODO

# Installation
include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# Macros from LibCEU.
set(LIBCEU_CMAKE_COLLECTION_DIR "${CMAKE_CURRENT_LIST_DIR}/deps/cmake_collections")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/enable_debug.cmake")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/libcmake/enhanced_find.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/boost.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/parallel.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/htslib.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/mpi.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/rng.cmake")
# OpenMP temporarily disabled. include("${CMAKE_CURRENT_LIST_DIR}/cmake_includes/openmp.cmake")

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/deps/slim_libceu")
include_directories("${CMAKE_CURRENT_LIST_DIR}/deps/slim_libceu")
include_directories("${CMAKE_CURRENT_LIST_DIR}/deps/concurrentqueue")
include_directories("${CMAKE_CURRENT_LIST_DIR}/deps/simde")

set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} libceu)
message(NOTICE "ART_MODERN_LINK_LIBS=${ART_MODERN_LINK_LIBS}")

# Main
if(CEU_CM_SHOULD_ENABLE_TEST)
    configure_file("${CMAKE_CURRENT_LIST_DIR}/src/test/test_adaptor.h.in" "${CMAKE_BINARY_DIR}/include/test_adaptor.h")
endif()
configure_file("${CMAKE_CURRENT_LIST_DIR}/src/lib/art_modern_config.h.in"
               "${CMAKE_BINARY_DIR}/include/art_modern_config.h")

include_directories("${CMAKE_CURRENT_LIST_DIR}/src/lib")
include_directories("${CMAKE_BINARY_DIR}/include/")
include_directories(${env.INCLUDE})
include_directories(${CMAKE_INCLUDE_PATH})

file(GLOB ART_MODERN_LIB_SOURCE "src/lib/*.cc" "src/lib/*/*.cc")
add_library(art_modern_lib ${ART_MODERN_LIB_SOURCE})

file(GLOB ART_SOURCE "src/art/*.cc")
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
    file(GLOB ART_MODERN_TEST_SOURCES "src/test/*.cc")

    foreach(ART_MODERN_TEST_SOURCE ${ART_MODERN_TEST_SOURCES})
        get_filename_component(ART_MODERN_TEST_EXEC "${ART_MODERN_TEST_SOURCE}" NAME_WLE)
        set(ART_MODERN_TEST_EXEC "am_test_${ART_MODERN_TEST_EXEC}")
        add_executable("${ART_MODERN_TEST_EXEC}" "${ART_MODERN_TEST_SOURCE}")
        target_link_libraries("${ART_MODERN_TEST_EXEC}" art_modern_lib Boost::unit_test_framework)
        add_test(NAME "${ART_MODERN_TEST_EXEC}" COMMAND "${ART_MODERN_TEST_EXEC}")
    endforeach()
else()
    message(NOTICE "Testing disabled since Boost::unit_test_framework was not found")
endif()

# if(CEU_CM_SHOULD_ENABLE_TEST) file(GLOB ART_MODERN_BENCHMARK_SOURCES "src/benchmark/*.cc")

# foreach(ART_MODERN_BENCHMARK_SOURCE ${ART_MODERN_BENCHMARK_SOURCES}) get_filename_component(ART_MODERN_BENCHMARK_EXEC
# "${ART_MODERN_BENCHMARK_SOURCE}" NAME_WLE) set(ART_MODERN_BENCHMARK_EXEC "am_bench_${ART_MODERN_BENCHMARK_EXEC}")
# add_executable("${ART_MODERN_BENCHMARK_EXEC}" "${ART_MODERN_BENCHMARK_SOURCE}")
# target_link_libraries("${ART_MODERN_BENCHMARK_EXEC}" art_modern_lib) endforeach() endif()

# file(GLOB ART_MODERN_UTIL_SOURCES "src/util/*.cc")

# foreach(ART_MODERN_UTIL_SOURCE ${ART_MODERN_UTIL_SOURCES}) get_filename_component(ART_MODERN_UTIL_EXEC
# "${ART_MODERN_UTIL_SOURCE}" NAME_WLE) set(ART_MODERN_UTIL_EXEC "am_util_${ART_MODERN_UTIL_EXEC}")
# add_executable("${ART_MODERN_UTIL_EXEC}" "${ART_MODERN_UTIL_SOURCE}") target_link_libraries("${ART_MODERN_UTIL_EXEC}"
# art_modern_lib) endforeach()
