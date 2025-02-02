cmake_minimum_required(VERSION 3.11)
project(
    test_helloworld_link
    VERSION 0.1
    LANGUAGES C)

include("${CMAKE_CURRENT_LIST_DIR}/../libcmake/set_static_target.cmake")

include_directories("${CMAKE_CURRENT_LIST_DIR}/../test_link_proj_source")
add_library(main_lib STATIC "${CMAKE_CURRENT_LIST_DIR}/../test_link_proj_source/lib.c")
add_executable(main_exe "${CMAKE_CURRENT_LIST_DIR}/../test_link_proj_source/main.c")

ceu_cm_set_static_target(main_exe)
ceu_cm_set_static_target(main_lib)
target_link_libraries(main_exe main_lib)

enable_testing()
add_test(NAME main_exe COMMAND main_exe)
