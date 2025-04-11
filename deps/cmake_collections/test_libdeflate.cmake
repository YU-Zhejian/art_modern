include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

if(NOT TARGET CEU_CM_EFL::libdeflate_shared)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libdeflate_shared PKGCONF_NAME libdeflate LINKER_FLAG deflate)
endif()
if(NOT TARGET CEU_CM_EFL::libdeflate_static)
    ceu_cm_enhanced_find_library(
        STATIC
        OUTPUT_VARIABLE
        libdeflate_static
        PKGCONF_NAME
        libdeflate
        LINKER_FLAG
        deflate)
endif()
if(TARGET CEU_CM_EFL::libdeflate_shared)
    ceu_cm_enhanced_try_run(
        VARNAME
        LIBDEFLATE
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libdeflate.c"
        LINK_LIBRARIES
        CEU_CM_EFL::libdeflate_shared
        DEPENDS
        C_HELLOWORLD)
        else()
        set(CEU_CM_HAVE_WORKING_LIBDEFLATE_RUN_SHARED 127 CACHE INTERNAL "libdeflate not found.")
        set(CEU_CM_HAVE_WORKING_LIBDEFLATE_COMPILE_SHARED OFF CACHE INTERNAL "libdeflate not found.")
endif()
if(TARGET CEU_CM_EFL::libdeflate_static)
    ceu_cm_enhanced_try_run(
        STATIC
        VARNAME
        LIBDEFLATE
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libdeflate.c"
        LINK_LIBRARIES
        CEU_CM_EFL::libdeflate_static
        DEPENDS
        C_HELLOWORLD)
        else()
        set(CEU_CM_HAVE_WORKING_LIBDEFLATE_RUN_STATIC 127 CACHE INTERNAL "libdeflate not found.")
        set(CEU_CM_HAVE_WORKING_LIBDEFLATE_COMPILE_STATIC OFF CACHE INTERNAL "libdeflate not found.")
endif()
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED" ON)
    ceu_cm_print_test_status("libdeflate" LIBDEFLATE)
endif()
