include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")

ceu_cm_enhanced_try_run(VARNAME C_HELLOWORLD SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_helloworld.c")
ceu_cm_enhanced_try_run(STATIC VARNAME C_HELLOWORLD SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_helloworld.c")

if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED" ON)
    ceu_cm_print_test_status("helloworld (c)" C_HELLOWORLD)
endif()
