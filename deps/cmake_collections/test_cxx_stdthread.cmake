include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_cxx_helloworld.cmake")

if(NOT DEFINED Threads_FOUND)
    find_package(Threads)
endif()

ceu_cm_enhanced_try_run(
    VARNAME
    CXX_STDTHREAD
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_stdthread.cc"
    LINK_LIBRARIES
    Threads::Threads
    "${CEU_CM_LIBM_SHARED}"
    DEPENDS
    CXX_HELLOWORLD)
ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    CXX_STDTHREAD
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_stdthread.cc"
    LINK_LIBRARIES
    Threads::Threads
    "${CEU_CM_LIBM_STATIC}"
    DEPENDS
    CXX_HELLOWORLD)
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")

    ceu_cm_print_test_status("stdthread (cxx)" CXX_STDTHREAD)
endif()
