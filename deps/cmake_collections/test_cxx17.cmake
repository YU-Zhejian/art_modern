include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_cxx14.cmake")

ceu_cm_enhanced_try_run(VARNAME CXX_CXX17 SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_cxx17.cc" DEPENDS CXX_CXX14)

ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    CXX_CXX17
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_cxx17.cc"
    DEPENDS
    CXX_CXX14)
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        ON
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("CXX17" CXX_CXX17)
endif()
