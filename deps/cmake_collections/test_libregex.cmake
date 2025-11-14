# Would firstly test whether C program runs without libregex.
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")

ceu_cm_enhanced_try_run(VARNAME C_NO_LIBREGEX SRC_PATH "${CMAKE_CURRENT_LIST_DIR}/src/test_libregex.c" DEPENDS
                        C_HELLOWORLD)
ceu_cm_enhanced_try_run(
    STATIC
    VARNAME
    C_NO_LIBREGEX
    SRC_PATH
    "${CMAKE_CURRENT_LIST_DIR}/src/test_libregex.c"
    DEPENDS
    C_HELLOWORLD)

if(NOT TARGET CEU_CM_EFL::libregex_shared)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libregex_shared PKGCONF_NAME regex)
endif()
if(NOT TARGET CEU_CM_EFL::libregex_static)
    ceu_cm_enhanced_find_library(STATIC OUTPUT_VARIABLE libregex_static PKGCONF_NAME regex)
endif()
# Some times, the system regex needs also libiconv.

if(NOT TARGET CEU_CM_EFL::libiconv_shared)
    ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libiconv_shared PKGCONF_NAME iconv)
endif()
if(NOT TARGET CEU_CM_EFL::libiconv_static)
    ceu_cm_enhanced_find_library(STATIC OUTPUT_VARIABLE libiconv_static PKGCONF_NAME iconv)
endif()
unset(CEU_CM_REGEX_REQUIRES_ICONV_STATIC)
unset(CEU_CM_REGEX_REQUIRES_ICONV_SHARED)

if(TARGET CEU_CM_EFL::libregex_shared)
    ceu_cm_enhanced_try_run(
        VARNAME
        C_WITH_LIBREGEX
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libregex.c"
        DEPENDS
        C_HELLOWORLD
        LINK_LIBRARIES
        CEU_CM_EFL::libregex_shared)
    if(NOT CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_RUN_SHARED EQUAL 0 AND TARGET CEU_CM_EFL::libiconv_shared)
        set(CEU_CM_REGEX_REQUIRES_ICONV_SHARED
            ON
            CACHE INTERNAL "libregex needs libiconv")
        # Try again with libiconv
        ceu_cm_enhanced_try_run(
            VARNAME
            C_WITH_LIBREGEX_ICONV
            SRC_PATH
            "${CMAKE_CURRENT_LIST_DIR}/src/test_libregex.c"
            DEPENDS
            C_HELLOWORLD
            LINK_LIBRARIES
            CEU_CM_EFL::libregex_shared
            CEU_CM_EFL::libiconv_shared)
    else()
        set(CEU_CM_REGEX_REQUIRES_ICONV_SHARED
            OFF
            CACHE INTERNAL "libregex does not need libiconv")
        set(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_ICONV_COMPILE_SHARED
            OFF
            CACHE BOOL "not tested")
        set(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_ICONV_RUN_SHARED
            255
            CACHE INTERNAL "not tested")
    endif()
else()
    set(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_COMPILE_SHARED
        OFF
        CACHE BOOL "libregex not found")
    set(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_RUN_SHARED
        255
        CACHE INTERNAL "libregex not found")
endif()
if(TARGET CEU_CM_EFL::libregex_static)
    ceu_cm_enhanced_try_run(
        STATIC
        VARNAME
        C_WITH_LIBREGEX
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_libregex.c"
        DEPENDS
        C_HELLOWORLD
        LINK_LIBRARIES
        CEU_CM_EFL::libregex_static)
    if(NOT CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_RUN_STATIC EQUAL 0 AND TARGET CEU_CM_EFL::libiconv_static)
        set(CEU_CM_REGEX_REQUIRES_ICONV_STATIC
            ON
            CACHE INTERNAL "libregex needs libiconv")
        # Try again with libiconv
        ceu_cm_enhanced_try_run(
            STATIC
            VARNAME
            C_WITH_LIBREGEX_ICONV
            SRC_PATH
            "${CMAKE_CURRENT_LIST_DIR}/src/test_libregex.c"
            DEPENDS
            C_HELLOWORLD
            LINK_LIBRARIES
            CEU_CM_EFL::libregex_static
            CEU_CM_EFL::libiconv_static)
    else()
        set(CEU_CM_REGEX_REQUIRES_ICONV_STATIC
            OFF
            CACHE INTERNAL "libregex does not need libiconv")
        set(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_ICONV_COMPILE_STATIC
            OFF
            CACHE BOOL "not tested")
        set(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_ICONV_RUN_STATIC
            255
            CACHE INTERNAL "not tested")
    endif()
else()
    set(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_COMPILE_STATIC
        OFF
        CACHE BOOL "libregex not found")
    set(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_RUN_STATIC
        255
        CACHE INTERNAL "libregex not found")
endif()
if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    ceu_cm_print_test_status("libregex: without -lregex (c)" C_NO_LIBREGEX)
    ceu_cm_print_test_status("libregex: with -lregex (c)" C_WITH_LIBREGEX)
    ceu_cm_print_test_status("libregex: with -lregex and -liconv (c)" C_WITH_LIBREGEX_ICONV)
endif()
if(CEU_CM_HAVE_WORKING_C_NO_LIBREGEX_RUN_SHARED EQUAL 0 AND NOT CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBREGEX_SHARED
        ""
        CACHE INTERNAL "regex functions work without libregex")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBREGEX_SHARED
        "CEU_CM_EFL::libregex_shared"
        CACHE INTERNAL "regex functions work with libregex")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_ICONV_RUN_SHARED EQUAL 0)
    set(CEU_CM_LIBREGEX_SHARED
        "CEU_CM_EFL::libregex_shared;CEU_CM_EFL::libiconv_shared"
        CACHE INTERNAL "regex functions work with libregex and iconv")
else()
    set(CEU_CM_LIBREGEX_SHARED
        "CEU_CM_LIBREGEX_SHARED-NOTFOUND"
        CACHE INTERNAL "regex functions not working")
endif()

if(CEU_CM_HAVE_WORKING_C_NO_LIBREGEX_RUN_STATIC EQUAL 0 AND NOT CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_RUN_STATIC EQUAL 0)
    set(CEU_CM_LIBREGEX_STATIC
        ""
        CACHE INTERNAL "regex functions work without libregex")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_RUN_STATIC EQUAL 0)
    set(CEU_CM_LIBREGEX_STATIC
        "CEU_CM_EFL::libregex_static"
        CACHE INTERNAL "regex functions work with libregex")
elseif(CEU_CM_HAVE_WORKING_C_WITH_LIBREGEX_ICONV_RUN_STATIC EQUAL 0)
    set(CEU_CM_LIBREGEX_STATIC
        "CEU_CM_EFL::libregex_static;CEU_CM_EFL::libiconv_static"
        CACHE INTERNAL "regex functions work with libregex and iconv")
else()
    set(CEU_CM_LIBREGEX_STATIC
        "CEU_CM_LIBREGEX_STATIC-NOTFOUND"
        CACHE INTERNAL "regex functions not working")
endif()
