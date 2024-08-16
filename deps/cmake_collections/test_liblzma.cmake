include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

if(NOT DEFINED LIBLZMA_LIBRARY_SHARED)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CEU_CM_LIBLZMA_LIBRARY_SHARED PKGCONFIG_NAME liblzma LINKER_FLAG
                                 lzma)
endif()

if(NOT DEFINED LIBLZMA_LIBRARY_STATIC)
    ceu_cm_enhanced_find_library(
        STATIC
        OUTPUT_VARIABLE
        CEU_CM_LIBLZMA_LIBRARY_STATIC
        PKGCONFIG_NAME
        liblzma
        LINKER_FLAG
        lzma)
endif()

ceu_cm_enhanced_try_run(
    VARNAME
    LIBLZMA
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_liblzma.c"
    LINK_LIBRARIES
    "${CEU_CM_LIBLZMA_LIBRARY_SHARED}"
    DEPENDS
    C_HELLOWORLD)
ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    LIBLZMA
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_liblzma.c"
    LINK_LIBRARIES
    "${CEU_CM_LIBLZMA_LIBRARY_STATIC}"
    DEPENDS
    C_HELLOWORLD)
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("liblzma" LIBLZMA)
endif()
