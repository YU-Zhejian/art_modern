include(FetchContent)
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_find.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/enhanced_try_run.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/print_test_status.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/test_cxx_stdthread.cmake")

if(NOT DEFINED CEU_CM_TRY_COMPILE_INCLUDES)
    set(CEU_CM_TRY_COMPILE_INCLUDES
        ""
        CACHE INTERNAL "Extra path added to try_compile")
endif()

function(ceu_cm_find_or_get_gtest)
    if(NOT DEFINED CEU_CM_FIND_OR_GET_GTEST_RESULT)
        set(CEU_CM_FIND_OR_GET_GTEST_RESULT
            FALSE
            CACHE INTERNAL "Default to fail.")
    endif()
    if(NOT CEU_CM_FIND_OR_GET_GTEST_RESULT)
        set(options "")
        set(oneValueArgs POLICY GIT_REPOSITORY GIT_TAG LOCAL_PREFIX)
        set(multiValueArgs "")
        cmake_parse_arguments(CEU_CM_FIND_OR_GET_GTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

        if(NOT DEFINED CEU_CM_FIND_OR_GET_GTEST_POLICY)
            set(CEU_CM_FIND_OR_GET_GTEST_POLICY "System")
        elseif("${CEU_CM_FIND_OR_GET_GTEST_POLICY}" STREQUAL "Download")
            if(NOT DEFINED CEU_CM_FIND_OR_GET_GTEST_GIT_REPOSITORY)
                set(CEU_CM_FIND_OR_GET_GTEST_GIT_REPOSITORY "https://github.com/google/googletest.git")
            endif()

            if(NOT DEFINED CEU_CM_FIND_OR_GET_GTEST_GIT_TAG)
                set(CEU_CM_FIND_OR_GET_GTEST_GIT_TAG "703bd9caab50b139428cea1aaff9974ebee5742e")
            endif()

            FetchContent_Declare(
                googletest
                GIT_REPOSITORY "${CEU_CM_FIND_OR_GET_GTEST_GIT_REPOSITORY}"
                GIT_TAG "${CEU_CM_FIND_OR_GET_GTEST_GIT_TAG}")
            if(NOT DEFINED googletest_POPULATED)
                message(
                    STATUS
                        "CEU_CM: Fetching GTest ${CEU_CM_FIND_OR_GET_GTEST_GIT_TAG}@${CEU_CM_FIND_OR_GET_GTEST_GIT_REPOSITORY}"
                )
                FetchContent_Populate(googletest)
            endif()
            set(CEU_CM_GTEST_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/usr)

            file(MAKE_DIRECTORY "${CEU_CM_GTEST_INSTALL_PREFIX}")
            file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/build_thirdparty")
            file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/build_thirdparty")

            if(NOT DEFINED CEU_CM_GTEST_CONFIGURE)
                if(CEU_CM_HAVE_WORKING_CXX_HELLOWORLD_RUN_SHARED EQUAL 0)
                    # TODO: Not implemented
                    message(
                        STATUS
                            "CEU_CM: Configuring SHARED GTest ${googletest_SOURCE_DIR} -> ${CMAKE_BINARY_DIR}/build_thirdparty"
                    )
                    execute_process(
                        COMMAND "cmake" -G "${CMAKE_GENERATOR}" -DCMAKE_BUILD_TYPE=Release "${googletest_SOURCE_DIR}"
                        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build_thirdparty"
                        ERROR_FILE "${CMAKE_BINARY_DIR}/compile_logs/gtest_configure.stderr.log"
                        OUTPUT_FILE "${CMAKE_BINARY_DIR}/compile_logs/gtest_configure.stdout.log"
                        RESULT_VARIABLE CEU_CM_GTEST_CONFIGURE)
                elseif(CEU_CM_HAVE_WORKING_CXX_HELLOWORLD_RUN_STATIC EQUAL 0)
                    message(
                        STATUS
                            "CEU_CM: Configuring STATIC GTest ${googletest_SOURCE_DIR} -> ${CMAKE_BINARY_DIR}/build_thirdparty"
                    )
                    execute_process(
                        COMMAND "cmake" -G "${CMAKE_GENERATOR}" -DCMAKE_BUILD_TYPE=Release "${googletest_SOURCE_DIR}"
                        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build_thirdparty"
                        ERROR_FILE "${CMAKE_BINARY_DIR}/compile_logs/gtest_configure.stderr.log"
                        OUTPUT_FILE "${CMAKE_BINARY_DIR}/compile_logs/gtest_configure.stdout.log"
                        RESULT_VARIABLE CEU_CM_GTEST_CONFIGURE)
                else()
                    message(
                        FATAL_ERROR
                            "CEU_CM: Building of GTest needs CEU_CM_HAVE_WORKING_CXX_HELLOWORLD_RUN_STATIC or CEU_CM_HAVE_WORKING_CXX_HELLOWORLD_RUN_SHARED"
                    )
                endif()

                if(NOT CEU_CM_GTEST_CONFIGURE EQUAL 0)
                    message(FATAL_ERROR "GTest configure failed -- ${CEU_CM_GTEST_CONFIGURE}!")
                endif()
            endif()

            if(NOT DEFINED CEU_CM_GTEST_BUILD)
                message(
                    STATUS
                        "CEU_CM: Building GTest ${googletest_SOURCE_DIR}/build_thirdparty -> ${CMAKE_BINARY_DIR}/build_thirdparty"
                )
                execute_process(
                    COMMAND "cmake" --build "${CMAKE_BINARY_DIR}/build_thirdparty"
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build_thirdparty"
                    ERROR_FILE "${CMAKE_BINARY_DIR}/compile_logs/gtest_build.stderr.log"
                    OUTPUT_FILE "${CMAKE_BINARY_DIR}/compile_logs/gtest_build.stdout.log"
                    RESULT_VARIABLE CEU_CM_GTEST_BUILD)

                if(NOT CEU_CM_GTEST_BUILD EQUAL 0)
                    message(FATAL_ERROR "CEU_CM: GTest build failed -- ${CEU_CM_GTEST_BUILD}!")
                endif()
            endif()

            if(NOT DEFINED CEU_CM_GTEST_INSTALL)
                message(
                    STATUS
                        "CEU_CM: Installing GTest ${googletest_SOURCE_DIR}/build_thirdparty -> ${CEU_CM_GTEST_INSTALL_PREFIX}"
                )
                execute_process(
                    COMMAND "cmake" --install "${CMAKE_BINARY_DIR}/build_thirdparty" --prefix
                            "${CEU_CM_GTEST_INSTALL_PREFIX}"
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build_thirdparty"
                    ERROR_FILE "${CMAKE_BINARY_DIR}/compile_logs/gtest_install.stderr.log"
                    OUTPUT_FILE "${CMAKE_BINARY_DIR}/compile_logs/gtest_install.stdout.log"
                    RESULT_VARIABLE CEU_CM_GTEST_INSTALL)

                if(NOT CEU_CM_GTEST_INSTALL EQUAL 0)
                    message(FATAL_ERROR "CEU_CM: GTest install failed -- ${CEU_CM_GTEST_INSTALL}!")
                endif()
            endif()

            list(APPEND CMAKE_PREFIX_PATH "${CEU_CM_GTEST_INSTALL_PREFIX}")
            list(APPEND CEU_CM_TRY_COMPILE_INCLUDES "${CEU_CM_GTEST_INSTALL_PREFIX}/include")
            include_directories(BEFORE "${CEU_CM_GTEST_INSTALL_PREFIX}/include")

        elseif("${CEU_CM_FIND_OR_GET_GTEST_POLICY}" STREQUAL "Local")
            if(NOT DEFINED CEU_CM_FIND_OR_GET_GTEST_LOCAL_PREFIX)
                message(FATAL_ERROR "CEU_CM: In local mode please define LOCAL_PREFIX!")
            endif()

            list(APPEND CMAKE_PREFIX_PATH "${CEU_CM_FIND_OR_GET_GTEST_LOCAL_PREFIX}")
            include_directories(BEFORE "${CEU_CM_FIND_OR_GET_GTEST_LOCAL_PREFIX}/include")
            list(APPEND CEU_CM_TRY_COMPILE_INCLUDES "${CEU_CM_FIND_OR_GET_GTEST_LOCAL_PREFIX}/include")
        endif()
        if(NOT MSVC)
            ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CEU_CM_GTEST_LIBRARY_SHARED PKGCONFIG_NAME gtest LINKER_FLAG
                                         gtest)
            ceu_cm_enhanced_find_library(
                STATIC
                OUTPUT_VARIABLE
                CEU_CM_GTEST_LIBRARY_STATIC
                PKGCONFIG_NAME
                gtest
                LINKER_FLAG
                gtest)
            ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CEU_CM_GTEST_MAIN_LIBRARY_SHARED PKGCONFIG_NAME gtest_main
                                         LINKER_FLAG gtest_main)
            ceu_cm_enhanced_find_library(
                STATIC
                OUTPUT_VARIABLE
                CEU_CM_GTEST_MAIN_LIBRARY_STATIC
                PKGCONFIG_NAME
                gtest_main
                LINKER_FLAG
                gtest_main)
        else()
            ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CEU_CM_GTEST_LIBRARY_SHARED LINKER_FLAG gtestg)
            ceu_cm_enhanced_find_library(STATIC OUTPUT_VARIABLE CEU_CM_GTEST_LIBRARY_STATIC LINKER_FLAG gtestd)
            ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CEU_CM_GTEST_MAIN_LIBRARY_SHARED LINKER_FLAG gtest_maind)
            ceu_cm_enhanced_find_library(STATIC OUTPUT_VARIABLE CEU_CM_GTEST_MAIN_LIBRARY_STATIC LINKER_FLAG
                                         gtest_maind)
        endif()

        set(CEU_CM_GTEST_MAIN_LIBRARY_SHARED ${CEU_CM_GTEST_MAIN_LIBRARY_SHARED} ${CEU_CM_GTEST_LIBRARY_SHARED}
                                             Threads::Threads)
        set(CEU_CM_GTEST_MAIN_LIBRARY_STATIC ${CEU_CM_GTEST_MAIN_LIBRARY_STATIC} ${CEU_CM_GTEST_LIBRARY_STATIC}
                                             Threads::Threads)

        ceu_cm_enhanced_try_run(
            VARNAME
            GTEST
            SRC_PATH
            "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/src/test_gtest.cc"
            LINK_LIBRARIES
            ${CEU_CM_GTEST_MAIN_LIBRARY_SHARED}
            DEPENDS
            CXX_STDTHREAD)
        ceu_cm_enhanced_try_run(
            STATIC
            VARNAME
            GTEST
            SRC_PATH
            "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/src/test_gtest.cc"
            LINK_LIBRARIES
            ${CEU_CM_GTEST_MAIN_LIBRARY_STATIC}
            DEPENDS
            CXX_STDTHREAD)

        ceu_cm_print_test_status("gtest (cxx)" GTEST "GTEST_MAIN")
        if(CEU_CM_HAVE_WORKING_GTEST_RUN_SHARED EQUAL 0 OR CEU_CM_HAVE_WORKING_GTEST_RUN_STATIC EQUAL 0)
            set(CEU_CM_FIND_OR_GET_GTEST_RESULT
                TRUE
                CACHE INTERNAL "GTest found.")
        endif()
    endif()
endfunction()
