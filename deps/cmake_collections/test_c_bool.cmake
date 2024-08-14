include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

ceu_cm_enhanced_try_run(
        VARNAME
        C_UNDERLINE_BOOL_KEYWORD
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_c_bool.c"
        DEPENDS
        C_HELLOWORLD
        COMPILE_DEFS
        C_UNDERLINE_BOOL_KEYWORD
)

ceu_cm_enhanced_try_run(
        STATIC
        VARNAME
        C_UNDERLINE_BOOL_KEYWORD
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_c_bool.c"
        DEPENDS
        C_HELLOWORLD
        COMPILE_DEFS
        C_UNDERLINE_BOOL_KEYWORD
)

ceu_cm_enhanced_try_run(
        VARNAME
        C_STDBOOL
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_c_bool.c"
        DEPENDS
        C_HELLOWORLD
        COMPILE_DEFS
        C_STDBOOL
)

ceu_cm_enhanced_try_run(
        STATIC
        VARNAME
        C_STDBOOL
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_c_bool.c"
        DEPENDS
        C_HELLOWORLD
        COMPILE_DEFS
        C_STDBOOL
)

if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("_Bool (c)" C_UNDERLINE_BOOL_KEYWORD)
    ceu_cm_print_test_status("<stdbool.h> + bool (c)" C_STDBOOL)
endif()
