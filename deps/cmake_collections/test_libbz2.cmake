include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

if(NOT TARGET CEU_CM_EFL::libbz2_shared)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libbz2_shared PKGCONF_NAME bzip2 LINKER_FLAG bz2)
endif()
if(NOT TARGET CEU_CM_EFL::libbz2_static)
    ceu_cm_enhanced_find_library(
        STATIC
        OUTPUT_VARIABLE
        libbz2_static
        PKGCONF_NAME
        bzip2
        LINKER_FLAG
        bz2)
endif()
if(TARGET CEU_CM_EFL::libbz2_shared)
    ceu_cm_enhanced_try_run(
        VARNAME
        LIBBZ2
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libbz2.c"
        LINK_LIBRARIES
        CEU_CM_EFL::libbz2_shared
        DEPENDS
        C_HELLOWORLD)
else()
    set(CEU_CM_HAVE_WORKING_LIBBZ2_RUN_SHARED
        127
        CACHE INTERNAL "libbz2 not found.")
    set(CEU_CM_HAVE_WORKING_LIBBZ2_COMPILE_SHARED
        OFF
        CACHE INTERNAL "libbz2 not found.")
endif()
if(TARGET CEU_CM_EFL::libbz2_static)
    ceu_cm_enhanced_try_run(
        STATIC
        VARNAME
        LIBBZ2
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libbz2.c"
        LINK_LIBRARIES
        CEU_CM_EFL::libbz2_static
        DEPENDS
        C_HELLOWORLD)
else()
    set(CEU_CM_HAVE_WORKING_LIBBZ2_RUN_STATIC
        127
        CACHE INTERNAL "libbz2 not found.")
    set(CEU_CM_HAVE_WORKING_LIBBZ2_COMPILE_STATIC
        OFF
        CACHE INTERNAL "libbz2 not found.")
endif()
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED" ON)
    ceu_cm_print_test_status("libbz2" LIBBZ2)
endif()
