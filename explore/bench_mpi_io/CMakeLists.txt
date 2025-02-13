cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 17)

project(bench_mpi_io)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(MPI REQUIRED COMPONENTS C)

set(LINK_LIBS MPI::MPI_C)

add_executable(bench_one_io_per_job bench_one_io_per_job.cc)
target_link_libraries(bench_one_io_per_job PUBLIC ${LINK_LIBS})

add_executable(bench_nompi bench_nompi.cc)
