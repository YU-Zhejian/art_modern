cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED BOOST_CONFIG_PROVIDED_BY_BOOST)
    set(BOOST_CONFIG_PROVIDED_BY_BOOST ON)
endif()
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
    if(BOOST_CONFIG_PROVIDED_BY_BOOST)
        cmake_policy(SET CMP0167 NEW)
    else()
        cmake_policy(SET CMP0167 OLD)
    endif()
endif()

project(bench_rand)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Boost
find_package(Boost REQUIRED COMPONENTS random)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

# Find MKL
find_package(MKL REQUIRED)
if(MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIRS})
endif()

# Find GSL
find_package(GSL REQUIRED)
if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
endif()

# Find Abseil
find_package(absl REQUIRED)

set(LINK_LIBS MKL::MKL GSL::gsl absl::random_random absl::random_distributions)

add_executable(bench_common_task bench_common_task.cc)
target_link_libraries(bench_common_task PUBLIC ${LINK_LIBS})

add_executable(bench_bits bench_bits.cc)
target_link_libraries(bench_bits PUBLIC ${LINK_LIBS})

add_executable(bench_dists bench_dists.cc)
target_link_libraries(bench_dists PUBLIC ${LINK_LIBS})

add_executable(bench_uniform_real bench_uniform_real.cc)
target_link_libraries(bench_uniform_real PUBLIC ${LINK_LIBS})
