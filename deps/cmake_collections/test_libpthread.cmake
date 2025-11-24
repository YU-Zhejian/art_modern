# Would firstly test whether C program runs without libpthread.
# Copyright 2024 Kitware, Inc. Distributed under the OSI-approved BSD 3-Clause License.
# From CMake 3.28
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_c_helloworld.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_libm.cmake")

add_library(CEU_CM_EFL::pthread_flag INTERFACE IMPORTED)
set_property(
    TARGET CEU_CM_EFL::pthread_flag
    PROPERTY INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:SHELL:-Xcompiler -pthread>"
             "$<$<AND:$<NOT:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>,$<NOT:$<COMPILE_LANGUAGE:Swift>>>:-pthread>"
             # Add -static flag if BUILD_SHARED_LIBS is OFF
)
set_property(
    TARGET CEU_CM_EFL::pthread_flag
    PROPERTY INTERFACE_LINK_OPTIONS "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:SHELL:-Xcompiler -pthread>"
             "$<$<AND:$<NOT:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>,$<NOT:$<COMPILE_LANGUAGE:Swift>>>:-pthread>"

    )

if(BUILD_SHARED_LIBS)
    ceu_cm_enhanced_try_run(
        VARNAME
        C_WITH_PTHREAD_FLAG
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_pthread.c"
        DEPENDS
        C_HELLOWORLD
        LINK_LIBRARIES
        CEU_CM_EFL::pthread_flag
        ${CEU_CM_LIBM_SHARED})
    if(CEU_CM_HAVE_WORKING_C_WITH_PTHREAD_FLAG_COMPILE_SHARED AND CEU_CM_HAVE_WORKING_C_WITH_PTHREAD_FLAG_RUN_SHARED
                                                                  EQUAL 0)
        set(CEU_CM_LIBPTHREAD_SHARED
            "CEU_CM_EFL::pthread_flag"
            CACHE BOOL "pthread flag works")
    endif()
else()
    ceu_cm_enhanced_try_run(
        VARNAME
        C_WITH_PTHREAD_FLAG
        SRC_PATH
        "${CMAKE_CURRENT_LIST_DIR}/src/test_pthread.c"
        DEPENDS
        STATIC
        C_HELLOWORLD
        LINK_LIBRARIES
        CEU_CM_EFL::pthread_flag
        ${CEU_CM_LIBM_STATIC})
    if(CEU_CM_HAVE_WORKING_C_WITH_PTHREAD_FLAG_COMPILE_STATIC AND CEU_CM_HAVE_WORKING_C_WITH_PTHREAD_FLAG_RUN_STATIC
                                                                  EQUAL 0)
        set(CEU_CM_LIBPTHREAD_STATIC
            "CEU_CM_EFL::pthread_flag"
            CACHE BOOL "pthread flag works")
    endif()
endif()

if(BUILD_SHARED_LIBS)
    if(NOT CEU_CM_HAVE_WORKING_C_WITH_PTHREAD_FLAG_RUN_SHARED EQUAL 0)
        if(NOT TARGET CEU_CM_EFL::libpthread_shared)
            ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libpthread_shared LINKER_FLAG pthread)
        endif()
        if(TARGET CEU_CM_EFL::libpthread_shared)
            ceu_cm_enhanced_try_run(
                VARNAME
                C_WITH_LIBPTHREAD
                SRC_PATH
                "${CMAKE_CURRENT_LIST_DIR}/src/test_pthread.c"
                DEPENDS
                C_HELLOWORLD
                LINK_LIBRARIES
                CEU_CM_EFL::libpthread_shared
                ${CEU_CM_LIBM_SHARED})
        else()
            set(CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_COMPILE_SHARED
                OFF
                CACHE BOOL "libpthread not found")
            set(CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_RUN_SHARED
                255
                CACHE INTERNAL "libpthread not found")
        endif()
        if(CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_COMPILE_SHARED AND CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_RUN_SHARED
                                                                    EQUAL 0)
            set(CEU_CM_LIBPTHREAD_SHARED
                "CEU_CM_EFL::libpthread_shared"
                CACHE BOOL "-lpthread works")
        else()
            set(CEU_CM_LIBPTHREAD_SHARED
                "CEU_CM_LIBPTHREAD_SHARED-NOTFOUND"
                CACHE BOOL "-lpthread does not work")
        endif()
    endif()
else()
    if(NOT CEU_CM_HAVE_WORKING_C_WITH_PTHREAD_FLAG_RUN_STATIC EQUAL 0)
        if(NOT TARGET CEU_CM_EFL::libpthread_static)
            ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libpthread_static LINKER_FLAG pthread)
        endif()
        if(TARGET CEU_CM_EFL::libpthread_static)
            ceu_cm_enhanced_try_run(
                VARNAME
                C_WITH_LIBPTHREAD
                SRC_PATH
                "${CMAKE_CURRENT_LIST_DIR}/src/test_pthread.c"
                DEPENDS
                C_HELLOWORLD
                STATIC
                LINK_LIBRARIES
                CEU_CM_EFL::libpthread_static
                ${CEU_CM_LIBM_STATIC})
        else()
            set(CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_COMPILE_STATIC
                OFF
                CACHE BOOL "libpthread not found")
            set(CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_RUN_STATIC
                255
                CACHE INTERNAL "libpthread not found")
        endif()
        if(CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_COMPILE_STATIC AND CEU_CM_HAVE_WORKING_C_WITH_LIBPTHREAD_RUN_STATIC
                                                                    EQUAL 0)
            set(CEU_CM_LIBPTHREAD_STATIC
                "CEU_CM_EFL::libpthread_static"
                CACHE BOOL "-lpthread works")
        else()
            set(CEU_CM_LIBPTHREAD_STATIC
                "CEU_CM_LIBPTHREAD_STATIC-NOTFOUND"
                CACHE BOOL "-lpthread does not work")
        endif()
    endif()
endif()

if(NOT DEFINED "${CMAKE_CURRENT_LIST_FILE}_INCLUDED")
    set("${CMAKE_CURRENT_LIST_FILE}_INCLUDED" ON)
    ceu_cm_print_test_status("libpthread (-lpthread)" C_WITH_LIBPTHREAD)
    ceu_cm_print_test_status("libpthread (-pthread)" C_WITH_PTHREAD_FLAG)
endif()
