include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c99.cmake")

ceu_cm_enhanced_try_run(VARNAME C_C11 SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_c11.c" DEPENDS C99)

ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C_C11
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_c11.c"
    DEPENDS
    C99)
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("C11" C_C11)
endif()
