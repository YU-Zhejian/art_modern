include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c89.cmake")

ceu_cm_enhanced_try_run(VARNAME C_C99 SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_c99.c" DEPENDS C_C89)

ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C_C99
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_c99.c"
    DEPENDS
    C_C89)
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("C99" C_C99)
endif()
