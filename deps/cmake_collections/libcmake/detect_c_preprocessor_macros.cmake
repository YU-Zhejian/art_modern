#[=======================================================================[
ceu_cm_detect_c_preprocessor_macros -- The wrapper of <shell/ceu_cm_detect_c_preprocessor_macros.sh>

This script wouild automatically check language enabling status by checking definition of
`CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`.

Synopsis: ceu_cm_detect_c_preprocessor_macros()

Requires:
    - CMake Variable `CMAKE_C_COMPILER`
    - CMake Variable `CMAKE_CXX_COMPILER`

FIXME: Failed in NVHPC
#]=======================================================================]
function(ceu_cm_detect_c_preprocessor_macros)
    if(WIN32)
        set(NULDEV NUL) # The Null device, which truncates all inpit
    else()
        set(NULDEV /dev/null)
    endif()
    if(DEFINED CMAKE_C_COMPILER)
        execute_process(COMMAND sh "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/shell/detect_c_preprocessor_macros.sh"
                                "${CMAKE_C_COMPILER}" c "${CMAKE_BINARY_DIR}/compile_logs/")
    endif()
    if(DEFINED CMAKE_CXX_COMPILER)
        execute_process(COMMAND sh "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/shell/detect_c_preprocessor_macros.sh"
                                "${CMAKE_CXX_COMPILER}" c++ "${CMAKE_BINARY_DIR}/compile_logs/")
    endif()
    unset(NULDEV)
endfunction()
