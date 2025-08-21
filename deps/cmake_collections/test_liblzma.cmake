include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

if(NOT TARGET CEU_CM_EFL::liblzma_shared)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE liblzma_shared PKGCONF_NAME liblzma LINKER_FLAG lzma)
endif()
if(NOT TARGET CEU_CM_EFL::liblzma_static)
    ceu_cm_enhanced_find_library(
        STATIC
        OUTPUT_VARIABLE
        liblzma_static
        PKGCONF_NAME
        liblzma
        LINKER_FLAG
        lzma)
endif()
if(TARGET CEU_CM_EFL::liblzma_shared)
    ceu_cm_enhanced_try_run(
        VARNAME
        LIBLZMA
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_liblzma.c"
        LINK_LIBRARIES
        CEU_CM_EFL::liblzma_shared
        DEPENDS
        C_HELLOWORLD)
else()
    set(CEU_CM_HAVE_WORKING_LIBLZMA_RUN_SHARED
        127
        CACHE INTERNAL "liblzma not found.")
    set(CEU_CM_HAVE_WORKING_LIBLZMA_COMPILE_SHARED
        OFF
        CACHE INTERNAL "liblzma not found.")
endif()
if(TARGET CEU_CM_EFL::liblzma_static)
    ceu_cm_enhanced_try_run(
        STATIC
        VARNAME
        LIBLZMA
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_liblzma.c"
        LINK_LIBRARIES
        CEU_CM_EFL::liblzma_static
        DEPENDS
        C_HELLOWORLD)
else()
    set(CEU_CM_HAVE_WORKING_LIBLZMA_RUN_STATIC
        127
        CACHE INTERNAL "liblzma not found.")
    set(CEU_CM_HAVE_WORKING_LIBLZMA_COMPILE_STATIC
        OFF
        CACHE INTERNAL "liblzma not found.")
endif()
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED" ON)
    ceu_cm_print_test_status("liblzma" LIBLZMA)
endif()
