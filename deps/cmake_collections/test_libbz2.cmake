include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

if(NOT DEFINED CMAKE_CM_LIBBZ2_LIBRARY_SHARED)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CMAKE_CM_LIBBZ2_LIBRARY_SHARED LINKER_FLAG bz2)
endif()

if(NOT DEFINED CMAKE_CM_LIBBZ2_LIBRARY_STATIC)
    ceu_cm_enhanced_find_library(STATIC OUTPUT_VARIABLE CMAKE_CM_LIBBZ2_LIBRARY_STATIC LINKER_FLAG bz2)
endif()

ceu_cm_enhanced_try_run(
    VARNAME
    LIBBZ2
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_libbz2.c"
    LINK_LIBRARIES
    "${CMAKE_CM_LIBBZ2_LIBRARY_SHARED}"
    DEPENDS
    C_HELLOWORLD)
ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    LIBBZ2
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_libbz2.c"
    LINK_LIBRARIES
    "${CMAKE_CM_LIBBZ2_LIBRARY_STATIC}"
    DEPENDS
    C_HELLOWORLD)
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED"
        TRUE
        CACHE INTERNAL "This file was included")
    ceu_cm_print_test_status("libbz2" LIBBZ2)
endif()
