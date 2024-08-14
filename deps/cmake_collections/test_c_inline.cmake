include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

ceu_cm_enhanced_try_run(VARNAME C___INLINE SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_c___inline.c" DEPENDS
                        C_HELLOWORLD)

ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C___INLINE
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_c___inline.c"
    DEPENDS
    C_HELLOWORLD)

ceu_cm_enhanced_try_run(VARNAME C_INLINE SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_c_inline.c" DEPENDS C_HELLOWORLD)

ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C_INLINE
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_c_inline.c"
    DEPENDS
    C_HELLOWORLD)

ceu_cm_enhanced_try_run(VARNAME C___INLINE__ SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_c___inline__.c" DEPENDS
                        C_HELLOWORLD)

ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C___INLINE__
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_c___inline__.c"
    DEPENDS
    C_HELLOWORLD)

if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("__inline (c)" C___INLINE)
    ceu_cm_print_test_status("__inline__ (c)" C___INLINE__)
    ceu_cm_print_test_status("inline (c)" C_INLINE)
endif()
