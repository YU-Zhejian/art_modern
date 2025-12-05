#[=======================================================================[
ceu_cm_detect_c_preprocessor_macros -- The wrapper of <shell/ceu_cm_detect_c_preprocessor_macros.sh>

This script wouild automatically check language enabling status by checking definition of
`CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`.

Synopsis: ceu_cm_detect_c_preprocessor_macros()

Requires:
    - CMake Variable `CMAKE_C_COMPILER`
    - CMake Variable `CMAKE_CXX_COMPILER`
#]=======================================================================]

set(CEU_CM_DCPPM_MODULE_BASE_DIR "${CMAKE_CURRENT_LIST_DIR}")

function(ceu_cm_detect_c_preprocessor_macros)
    if(DEFINED CMAKE_C_COMPILER
       AND CMAKE_C_COMPILER
       AND CMAKE_C_COMPILER_ID)
        execute_process(COMMAND sh "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/shell/detect_c_preprocessor_macros.sh"
                                "${CMAKE_C_COMPILER}" c "${CMAKE_BINARY_DIR}/compile_logs/")
        execute_process(COMMAND sh "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/shell/detect_c_compiler_log.sh"
                                "${CMAKE_C_COMPILER}" c "${CMAKE_BINARY_DIR}/compile_logs/")
    endif()
    if(DEFINED CMAKE_CXX_COMPILER
       AND CMAKE_CXX_COMPILER
       AND CMAKE_CXX_COMPILER_ID)
        execute_process(COMMAND sh "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/shell/detect_c_preprocessor_macros.sh"
                                "${CMAKE_CXX_COMPILER}" c++ "${CMAKE_BINARY_DIR}/compile_logs/")
        execute_process(COMMAND sh "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/shell/detect_c_compiler_log.sh"
                                "${CMAKE_CXX_COMPILER}" c++ "${CMAKE_BINARY_DIR}/compile_logs/")
    endif()
endfunction()
