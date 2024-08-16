# Search for log4cpp <https://log4cpp.sourceforge.net>

include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_cxx_stdthread.cmake")

if(NOT DEFINED CEU_CM_LOG4CPP_LIBRARY_SHARED)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CEU_CM_LOG4CPP_LIBRARY_SHARED PKGCONFIG_NAME log4cpp LINKER_FLAG
                                 log4cpp)
endif()

if(NOT DEFINED CEU_CM_LOG4CPP_LIBRARY_STATIC)
    ceu_cm_enhanced_find_library(
        STATIC
        OUTPUT_VARIABLE
        CEU_CM_LOG4CPP_LIBRARY_STATIC
        PKGCONFIG_NAME
        log4cpp
        LINKER_FLAG
        log4cpp)
endif()

ceu_cm_enhanced_try_run(
    VARNAME
    LOG4CPP
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_log4cpp.cc"
    LINK_LIBRARIES
    ${CEU_CM_LOG4CPP_LIBRARY_SHARED}
    Threads::Threads
    DEPENDS
    CXX_STDTHREAD)

ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    LOG4CPP
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_log4cpp.cc"
    LINK_LIBRARIES
    ${CEU_CM_LOG4CPP_LIBRARY_SATTIC}
    Threads::Threads
    DEPENDS
    CXX_STDTHREAD)
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("log4cpp" LOG4CPP)
endif()
