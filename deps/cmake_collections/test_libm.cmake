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

if(NOT TARGET CEU_CM_EFL::libm_shared)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libm_shared LINKER_FLAG m)
endif()
if(NOT TARGET CEU_CM_EFL::libm_static)
    ceu_cm_enhanced_find_library(STATIC OUTPUT_VARIABLE libm_static LINKER_FLAG m)
endif()
if(TARGET CEU_CM_EFL::libm_shared)
    ceu_cm_enhanced_try_run(
        VARNAME
        C_WITH_LIBM
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libm.c"
        DEPENDS
        C_HELLOWORLD
        LINK_LIBRARIES
        CEU_CM_EFL::libm_shared)
else()

    set(CEU_CM_HAVE_WORKING_C_WITH_LIBM_COMPILE_SHARED
        OFF
        CACHE BOOL "libm not found")
    set(CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_SHARED
        255
        CACHE INTERNAL "libm not found")
endif()
if(TARGET CEU_CM_EFL::libm_static)
    ceu_cm_enhanced_try_run(
        STATIC
        VARNAME
        C_WITH_LIBM
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libm.c"
        DEPENDS
        C_HELLOWORLD
        LINK_LIBRARIES
        CEU_CM_EFL::libm_static)
else()
    set(CEU_CM_HAVE_WORKING_C_WITH_LIBM_COMPILE_STATIC
        OFF
        CACHE BOOL "libm not found")
    set(CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_STATIC
        255
        CACHE INTERNAL "libm not found")
endif()
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    ceu_cm_print_test_status("libm: without -lm (c)" C_NO_LIBM)
    ceu_cm_print_test_status("libm: with -lm (c)" C_WITH_LIBM)
endif()
if(CEU_CM_HAVE_WORKING_C_NO_LIBM_RUN_SHARED EQUAL 0 AND NOT CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBM_SHARED
        ""
        CACHE INTERNAL "mathematical functions work without libm")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBM_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBM_SHARED
        "CEU_CM_EFL::libm_shared"
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
        "CEU_CM_EFL::libm_static"
        CACHE INTERNAL "mathematical functions work with libm")
else()
    set(CEU_CM_LIBM_STATIC
        "CEU_CM_LIBM_STATIC-NOTFOUND"
        CACHE INTERNAL "mathematical functions not working")
endif()
