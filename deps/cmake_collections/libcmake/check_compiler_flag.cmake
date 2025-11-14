include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)

#[=======================================================================[
ceu_cm_enhanced_check_compiler_flag -- Check whether one of the flags in a group of flags is available.

Synopsis: ceu_cm_enhanced_check_compiler_flag(FLAG [[FLAG]...])

Params:
    - `FLAG`: A compiler flag.

Notice:
    - This function detects used languages by detecting `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER` CMake variable.
    - The options passed to this function is orderd. i.e., it would check whether the first flag is available. If so, the first flag would be added to OUT_VAR. If not, check the second flag until no flags left.
#]=======================================================================]
function(ceu_cm_enhanced_check_compiler_flag)
    set(options REQUIRED)
    set(oneValueArgs OUT_NAME)
    set(multiValueArgs FLAGS)
    cmake_parse_arguments(CEU_CM_ECCF "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    foreach(FLAG ${CEU_CM_ECCF_FLAGS})
        if(DEFINED CMAKE_C_COMPILER)
            if(NOT DEFINED C_COMPILER_HAVE_${FLAG})
                check_c_compiler_flag(${FLAG} C_COMPILER_HAVE_${FLAG})
                set(C_COMPILER_HAVE_${FLAG}
                    ${C_COMPILER_HAVE_${FLAG}}
                    CACHE INTERNAL "Whether C compiler supports ${FLAG}")
            endif()
        else()
            # If no C support is added, bypass the test.
            set(C_COMPILER_HAVE_${FLAG}
                ON
                CACHE INTERNAL "Test bypassed because no C compiler is defined.")
        endif()
        if(DEFINED CMAKE_CXX_COMPILER AND CMAKE_CXX_COMPILER)
            if(NOT DEFINED CXX_COMPILER_HAVE_${FLAG})
                check_cxx_compiler_flag(${FLAG} CXX_COMPILER_HAVE_${FLAG})
                set(CXX_COMPILER_HAVE_${FLAG}
                    ${CXX_COMPILER_HAVE_${FLAG}}
                    CACHE INTERNAL "Whether CXX compiler supports ${FLAG}")
            endif()
        else()
            set(CXX_COMPILER_HAVE_${FLAG}
                ON
                CACHE INTERNAL "Test bypassed because no CXX compiler is defined.")
        endif()
        if(C_COMPILER_HAVE_${FLAG} AND CXX_COMPILER_HAVE_${FLAG})
            # add_compile_options(${FLAG})
            set(${CEU_CM_ECCF_OUT_NAME}
                ${FLAG} ${${CEU_CM_ECCF_OUT_NAME}}
                PARENT_SCOPE)
            return()
        endif()
        unset(C_COMPILER_HAVE_${FLAG})
        unset(CXX_COMPILER_HAVE_${FLAG})
    endforeach()
    if(CEU_CM_ECCF_REQUIRED)
        message(FATAL_ERROR "None of the compiler flags ${CEU_CM_ECCF_FLAGS} are supported by the compiler.")
    endif()
endfunction()
