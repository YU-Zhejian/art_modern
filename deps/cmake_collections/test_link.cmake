if(NOT DEFINED CEU_CM_TRY_COMP_LINK_STATIC)
    execute_process(
        COMMAND
            ctest --verbose --stop-on-failure --build-and-test "${CMAKE_CURRENT_LIST_DIR}/test_link_proj_static"
            "${CMAKE_BINARY_DIR}/test_link_proj_static" --build-generator "${CMAKE_GENERATOR}" --test-command ctest -C
            Debug --build-options -DCMAKE_C_COMPILER="${CMAKE_C_COMPILER}" -DCMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
        RESULT_VARIABLE CEU_CM_TRY_COMP_LINK_STATIC_RETV
        OUTPUT_FILE "${CMAKE_BINARY_DIR}/compile_logs/test_link_proj_static_compile.log"
        ERROR_FILE "${CMAKE_BINARY_DIR}/compile_logs/test_link_proj_static_compile.err")

    if(NOT CEU_CM_TRY_COMP_LINK_STATIC_RETV EQUAL 0)
        set(CEU_CM_TRY_COMP_LINK_STATIC
            FALSE
            CACHE INTERNAL "doc")
    else()
        set(CEU_CM_TRY_COMP_LINK_STATIC
            TRUE
            CACHE INTERNAL "doc")
    endif()
endif()
if(NOT DEFINED CEU_CM_TRY_COMP_LINK_SHARED)
    execute_process(
        COMMAND
            ctest --verbose --stop-on-failure --build-and-test "${CMAKE_CURRENT_LIST_DIR}/test_link_proj_shared"
            "${CMAKE_BINARY_DIR}/test_link_proj_shared" --build-generator "${CMAKE_GENERATOR}" --test-command ctest -C
            Debug -DCMAKE_C_COMPILER="${CMAKE_C_COMPILER}" -DCMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
        RESULT_VARIABLE CEU_CM_TRY_COMP_LINK_SHARED_RETV
        OUTPUT_FILE "${CMAKE_BINARY_DIR}/compile_logs/test_link_proj_shared_compile.log"
        ERROR_FILE "${CMAKE_BINARY_DIR}/compile_logs/test_link_proj_shared_compile.err")
    if(NOT CEU_CM_TRY_COMP_LINK_SHARED_RETV EQUAL 0)
        set(CEU_CM_TRY_COMP_LINK_SHARED
            FALSE
            CACHE INTERNAL "doc")
    else()
        set(CEU_CM_TRY_COMP_LINK_SHARED
            TRUE
            CACHE INTERNAL "doc")
    endif()
endif()
