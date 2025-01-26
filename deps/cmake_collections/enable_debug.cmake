#[=======================================================================[
enable_debug.cmake -- C/C++ set up utilities

Requires:
    - `CMAKE_BUILD_TYPE`
    - `CEU_CM_SHOULD_ENABLE_TEST`
    - `CEU_CM_SHOULD_USE_NATIVE`

Function:
    - Defines several functions.
    - Check whether the CMake variable `CEU_CM_ENABLE_DEBUG_CMAKE_WAS_ALREADY_INCLUDED` was set. If so, skip all processes described below. if not, set this variable.
    - Detect C/C++ preprocessor macros. If this is not a C/C++ project, this step would be omitted.
    - Detect whether test should be built. If CMake variable `CEU_CM_SHOULD_ENABLE_TEST` was set, would skip this step. Otherwise, will setup test if the CMake variable `CMAKE_BUILD_TYPE` is not `Release`.
    - Detect whether we should build the application with native hardware support. If CMake variable `CEU_CM_SHOULD_USE_NATIVE` was not set, would set it to `OFF`. If this variable was set `ON`, would try and set `-march=native` `-mtune=native` `-mtune` flags.
    - Detect build type using CMake variable `CMAKE_BUILD_TYPE`.
        - If `Release`, will supress warnings (`-W0` `-w`), supress debug information (`-g0`) and enable optimization (`-Ofast` `-O3` `-O2`)
        - If `RelWithDebInfo` (Release with Debug Information), will supress warnings, enable debug information (`-g`) and enable optimization.
        - If `Debug` or if is others, will enable all warnings (`-Wall`) and extra warnings (`-Wextra`), enable pedantic mode (`-pedantic` `-Wpedantic`), enable full debug information (`-Og`) (`-g3`) and supress optimization (`-Og`). It will add `CEU_CM_IS_DEBUG` to compiler flags (i.e., `-DCEU_CM_IS_DEBUG`).
    - Print basic C/C++/Fortran information, including test and native information.
#]=======================================================================]

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/detect_c_preprocessor_macros.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/set_static_target.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/libcmake/check_compiler_flag.cmake")

#[=======================================================================[
ceu_cm_global_enhanced_check_compiler_flag -- Check whether one of the flags in a group of flags is available.

Synopsis: ceu_cm_global_enhanced_check_compiler_flag(FLAG [[FLAG]...])

Params:
    - `FLAG`: A linker flag.

Notice:
    - This function detects used languages by detecting `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER` CMake variable.
    - The options passed to this function is orderd. i.e., it would check whether the first flag is available. If so, the first flag would be added to compile options using CMake builtin function `add_compile_options`. If not, check the second flag until no flags left.

Sets:
    - CMake variable `C_COMPILER_HAVE_${FLAG}`.
#]=======================================================================]

if(NOT DEFINED CEU_CM_ENABLE_DEBUG_CMAKE_WAS_ALREADY_INCLUDED)
    # For debugging purposes, export compile commands.
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    set(CEU_CM_ENABLE_DEBUG_CMAKE_WAS_ALREADY_INCLUDED
        ON
        CACHE INTERNAL "Whether a description on environment was printed.")
    # Detect C/CXX Pre-Processor Macros
    if(NOT MSVC)
        ceu_cm_detect_c_preprocessor_macros()
    endif()

    # Detect build type.
    if(NOT DEFINED CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Debug")
    endif()

    # Detect Test.
    if(NOT DEFINED CEU_CM_SHOULD_ENABLE_TEST)
        if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
            set(CEU_CM_SHOULD_ENABLE_TEST
                ON
                CACHE INTERNAL "Test automatically enabled")
        else()
            set(CEU_CM_SHOULD_ENABLE_TEST
                OFF
                CACHE INTERNAL "Test automatically disabled")
        endif()
    endif()
    if(CEU_CM_SHOULD_ENABLE_TEST)
        enable_testing()
    endif()
    # Detect native.
    if(NOT DEFINED CEU_CM_SHOULD_USE_NATIVE)
        set(CEU_CM_SHOULD_USE_NATIVE OFF)
    endif()
endif()

if(NOT DEFINED ENV{CEU_CM_DEBUG_BANNER_SHOWN})
    set(ENV{CEU_CM_DEBUG_BANNER_SHOWN} 1)

    message(STATUS "/------------------- Basic Information -------------------\\")
    message(STATUS "|CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}, CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}")
    message(STATUS "|CMAKE_SYSTEM=${CMAKE_SYSTEM}")
    message(STATUS "|CMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}")
    message(STATUS "|CMAKE_AR=${CMAKE_AR}")
    message(STATUS "|CMAKE_RANLIB=${CMAKE_RANLIB}")
    message(STATUS "|CMAKE_EXECUTABLE_SUFFIX=${CMAKE_EXECUTABLE_SUFFIX}")
    message(
        STATUS
            "|CMAKE_SHARED_LIBRARY_PREFIX=${CMAKE_SHARED_LIBRARY_PREFIX}, CMAKE_SHARED_LIBRARY_SUFFIX=${CMAKE_SHARED_LIBRARY_SUFFIX}"
    )
    message(
        STATUS
            "|CMAKE_STATIC_LIBRARY_PREFIX=${CMAKE_STATIC_LIBRARY_PREFIX}, CMAKE_STATIC_LIBRARY_SUFFIX=${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )

    if(DEFINED CMAKE_C_COMPILER)
        message(
            STATUS
                "|CMAKE_C_COMPILER=${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID}|${CMAKE_C_COMPILER_ABI}) ver. ${CMAKE_C_COMPILER_VERSION} std. ${CMAKE_C_STANDARD}"
        )
        message(STATUS "|CMAKE_C_LINK_EXECUTABLE=${CMAKE_C_LINK_EXECUTABLE}")
    endif()
    if(DEFINED CMAKE_CXX_COMPILER)
        message(
            STATUS
                "|CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID}|${CMAKE_CXX_COMPILER_ABI}) ver. ${CMAKE_CXX_COMPILER_VERSION} std. ${CMAKE_CXX_STANDARD}"
        )
        message(STATUS "|CMAKE_CXX_LINK_EXECUTABLE=${CMAKE_CXX_LINK_EXECUTABLE}")
    endif()

    message(STATUS "|CMAKE_SYSTEM_LIBRARY_PATH=${CMAKE_SYSTEM_LIBRARY_PATH}")
    message(STATUS "|CMAKE_SYSTEM_INCLUDE_PATH=${CMAKE_SYSTEM_INCLUDE_PATH}")
    message(STATUS "|CMAKE_SYSTEM_PREFIX_PATH=${CMAKE_SYSTEM_PREFIX_PATH}")
    message(STATUS "|CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")
    message(STATUS "|CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

    message(
        STATUS
            "|CEU_CM_SHOULD_ENABLE_TEST=${CEU_CM_SHOULD_ENABLE_TEST}; CEU_CM_SHOULD_USE_NATIVE=${CEU_CM_SHOULD_USE_NATIVE}"
    )

    message(STATUS "\\------------------- Basic Information -------------------/")
endif()
