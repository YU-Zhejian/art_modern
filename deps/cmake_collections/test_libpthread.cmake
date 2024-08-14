# Would firstly test whether C program runs without libm.
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_libm.cmake")

ceu_cm_enhanced_try_run(
        VARNAME
        C_NO_LIBPTHREAD
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_pthread.c"
        DEPENDS
        C_HELLOWORLD
        LINK_LIBRARIES "${LIBM_LIBRARY_SHARED}"
)
ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
        C_NO_LIBPTHREAD
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_pthread.c"
    DEPENDS
    C_HELLOWORLD
        LINK_LIBRARIES "${LIBM_LIBRARY_STATIC}"
)

if(NOT DEFINED LIBPTHREAD_LIBRARY_SHARED)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE LIBPTHREAD_LIBRARY_SHARED LINKER_FLAG pthread)
endif()

if(NOT DEFINED LIBPTHREAD_LIBRARY_STATIC)
    ceu_cm_enhanced_find_library(STATIC OUTPUT_VARIABLE LIBPTHREAD_LIBRARY_STATIC LINKER_FLAG pthread)
endif()

ceu_cm_enhanced_try_run(
    VARNAME
    C_WITH_LIBPTHREAD
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_pthread.c"
    DEPENDS
    C_HELLOWORLD
    LINK_LIBRARIES
    "${LIBPTHREAD_LIBRARY_SHARED}"
        "${LIBM_LIBRARY_SHARED}"
)
ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C_WITH_LIBPTHREAD
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_pthread.c"
    DEPENDS
    C_HELLOWORLD
    LINK_LIBRARIES
    "${LIBPTHREAD_LIBRARY_STATIC}"
        "${LIBM_LIBRARY_STATIC}")
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("libpthread: without -lpthread (c)" C_NO_LIBPTHREAD)
    ceu_cm_print_test_status("libpthread: with -lpthread (c)" C_WITH_LIBPTHREAD)
endif()
if(CEU_CM_HAVE_WORKING_C_NO_LIBPTHREAD_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBPTHREAD_SHARED
        ""
        CACHE INTERNAL "pthread functions work without libm")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBPTHREAD_SHARED
        "${LIBPTHREAD_LIBRARY_SHARED}"
        CACHE INTERNAL "pthread functions work with libm")
else()
    set(CEU_CM_LIBPTHREAD_SHARED
        "CEU_CM_LIBPTHREAD_SHARED-NOTFOUND"
        CACHE INTERNAL "pthread functions not working")
endif()

if(CEU_CM_HAVE_WORKING_C_NO_LIBPTHREAD_RUN_STATIC EQUAL 0)
    set(CEU_CM_LIBPTHREAD_STATIC
        ""
        CACHE INTERNAL "pthread functions work without libm")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_RUN_STATIC EQUAL 0)
    set(CEU_CM_LIBPTHREAD_STATIC
        "${LIBPTHREAD_LIBRARY_STATIC}"
        CACHE INTERNAL "pthread functions work with libm")
else()
    set(CEU_CM_LIBPTHREAD_STATIC
        "CEU_CM_LIBPTHREAD_STATIC-NOTFOUND"
        CACHE INTERNAL "pthread functions not working")
endif()
