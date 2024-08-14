cmake_minimum_required(VERSION 3.17)
project(
    CPPTETGS_DIST_CMAKE
    VERSION 0.2.1
    LANGUAGES C CXX Fortran)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_FORTRAN_STANDARD 90)
set(CMAKE_WARN_DEPRECATED TRUE)

file(GLOB CPPTETGS_DIST_CMAKE_TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/test_*.cmake")

foreach(fiename ${CPPTETGS_DIST_CMAKE_TEST_FILES})
    include(${fiename})
endforeach()

file(GLOB CPPTETGS_DIST_CMAKE_FIND_OR_GET_FILES "${CMAKE_CURRENT_LIST_DIR}/find_or_get_*.cmake")

foreach(fiename ${CPPTETGS_DIST_CMAKE_FIND_OR_GET_FILES})
    include(${fiename})
endforeach()
