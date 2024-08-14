include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

ceu_cm_enhanced_try_run(VARNAME C___ATTRIBUTE__ SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_c___attribute__.c" DEPENDS
                        C_HELLOWORLD)

ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C___ATTRIBUTE__
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_c___attribute__.c"
    DEPENDS
    C_HELLOWORLD)

if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("__attribute__ (c)" C___ATTRIBUTE__)
endif()
