# Would firstly test whether C program runs without libm.
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

ceu_cm_enhanced_try_run(VARNAME C_NO_LIBM SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_libm.c" DEPENDS C_HELLOWORLD)
ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C_NO_LIBM
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_libm.c"
    DEPENDS
    C_HELLOWORLD)

if(NOT DEFINED LIBM_LIBRARY_SHARED)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE LIBM_LIBRARY_SHARED LINKER_FLAG m)
endif()

if(NOT DEFINED LIBM_LIBRARY_STATIC)
    ceu_cm_enhanced_find_library(STATIC OUTPUT_VARIABLE LIBM_LIBRARY_STATIC LINKER_FLAG m)
endif()

ceu_cm_enhanced_try_run(
    VARNAME
    C_WITH_LIBM
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_libm.c"
    DEPENDS
    C_HELLOWORLD
    LINK_LIBRARIES
    ${LIBM_LIBRARY_SHARED})
ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C_WITH_LIBM
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_libm.c"
    DEPENDS
    C_HELLOWORLD
    LINK_LIBRARIES
    ${LIBM_LIBRARY_STATIC})
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("libm: without -lm (c)" C_NO_LIBM)
    ceu_cm_print_test_status("libm: with -lm (c)" C_WITH_LIBM)
endif()
if(CEU_CM_HAVE_WORKING_C_NO_LIBM_RUN_SHARED EQUAL 0 AND NOT CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBM_SHARED
        ""
        CACHE INTERNAL "mathematical functions work without libm")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBM_SHARED
        "${LIBM_LIBRARY_SHARED}"
        CACHE INTERNAL "mathematical functions work with libm")
else()
    set(CEU_CM_LIBM_SHARED
        "CEU_CM_LIBM_SHARED-NOTFOUND"
        CACHE INTERNAL "mathematical functions not working")
endif()

if(CEU_CM_HAVE_WORKING_C_NO_LIBM_RUN_STATIC EQUAL 0 AND NOT CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_STATIC EQUAL 0)
    set(CEU_CM_LIBM_STATIC
        ""
        CACHE INTERNAL "mathematical functions work without libm")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_STATIC EQUAL 0)
    set(CEU_CM_LIBM_STATIC
        "${LIBM_LIBRARY_STATIC}"
        CACHE INTERNAL "mathematical functions work with libm")
else()
    set(CEU_CM_LIBM_STATIC
        "CEU_CM_LIBM_STATIC-NOTFOUND"
        CACHE INTERNAL "mathematical functions not working")
endif()
