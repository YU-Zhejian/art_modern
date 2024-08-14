cmake_minimum_required(VERSION 3.11)
project(
    test_helloworld_link
    VERSION 0.1
    LANGUAGES C)

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
add_compile_options(${CEU_CM_ADDITIONAL_COMPILER_FLAGS})

include_directories("${CMAKE_CURRENT_LIST_DIR}/../test_link_proj_source")
add_library(main_lib STATIC "${CMAKE_CURRENT_LIST_DIR}/../test_link_proj_source/lib.c")
add_executable(main_exe "${CMAKE_CURRENT_LIST_DIR}/../test_link_proj_source/main.c")

ceu_cm_set_static_target(main_exe)
ceu_cm_set_static_target(main_lib)
target_link_libraries(main_exe main_lib)

enable_testing()
add_test(NAME main_exe COMMAND main_exe)
