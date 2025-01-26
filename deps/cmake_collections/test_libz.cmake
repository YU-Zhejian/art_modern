# Search for zlib <https://www.zlib.net>

include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

if(NOT TARGET CEU_CM_EFL::libz_shared)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libz_shared PKGCONF_NAME zlib LINKER_FLAG z)
endif()
if(NOT TARGET CEU_CM_EFL::libz_static)
    ceu_cm_enhanced_find_library(
        STATIC
        OUTPUT_VARIABLE
        libz_static
        PKGCONF_NAME
        zlib
        LINKER_FLAG
        z)
endif()
if(TARGET CEU_CM_EFL::libz_shared)
    ceu_cm_enhanced_try_run(
        VARNAME
        LIBZ
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libz.c"
        LINK_LIBRARIES
        CEU_CM_EFL::libz_shared
        DEPENDS
        C_HELLOWORLD)
endif()
if(TARGET CEU_CM_EFL::libz_static)
    ceu_cm_enhanced_try_run(
        STATIC
        VARNAME
        LIBZ
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libz.c"
        LINK_LIBRARIES
        CEU_CM_EFL::libz_static
        DEPENDS
        C_HELLOWORLD)
endif()
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED" ON)
    ceu_cm_print_test_status("libz" LIBZ)
endif()
