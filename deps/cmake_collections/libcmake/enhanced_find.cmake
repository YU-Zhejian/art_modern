# Finds packages.

find_package(PkgConfig)

#[=======================================================================[
ceu_cm_get_abspath_from_linker_flag -- Get absolute path of libraries from linker flag.

Synopsis: `ceu_cm_get_abspath_from_linker_flag(OUTPUT_VARIABLE LINKER_FLAG IS_STATIC)`

Params:
    - `OUTPUT_VARIABLE`: Name of the output variable.
    - `LINKER_FLAG`: Linker flag without prefixes like `-l`.
    - `IS_STATIC`: Whether to find static or dynamic libraries.

Sample:
    - `ceu_cm_get_abspath_from_linker_flag(OV z TRUE)` -> `libz.a`
    - `ceu_cm_get_abspath_from_linker_flag(OV z FALSE)` -> `libz.so`

Sets:
    - `OUTPUT_VARIABLE`: Parent scope level.
    - `CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${IS_STATIC}`: Cache level. For debug purposes.

Warnings:
    - Finding libraries using linker flag cannot guarantee usability!

        For example, we have `app` -> `liba.so` -> `libb.so`.

        If `liba.so` does not contain implementations defined in `libb.so`,
        direct execution of `app` would yield bugs.
#]=======================================================================]
function(ceu_cm_get_abspath_from_linker_flag OUTPUT_VARIABLE LINKER_FLAG IS_STATIC)
    if(IS_STATIC)
        set(POSSIBLE_LIBRARY_FILENAMES
            "${LINKER_FLAG}${CMAKE_STATIC_LIBRARY_SUFFIX}"
            "${CMAKE_STATIC_LIBRARY_PREFIX}${LINKER_FLAG}${CMAKE_STATIC_LIBRARY_SUFFIX}"
            "lib${LINKER_FLAG}${CMAKE_STATIC_LIBRARY_SUFFIX}"
            "${LINKER_FLAG}.a"
            "${CMAKE_STATIC_LIBRARY_PREFIX}${LINKER_FLAG}.a"
            "lib${LINKER_FLAG}.a")
    else()
        set(POSSIBLE_LIBRARY_FILENAMES
            "${LINKER_FLAG}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${CMAKE_SHARED_LIBRARY_PREFIX}${LINKER_FLAG}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "lib${LINKER_FLAG}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${LINKER_FLAG}.so"
            "${CMAKE_SHARED_LIBRARY_PREFIX}${LINKER_FLAG}.so"
            "lib${LINKER_FLAG}.so"
            "${LINKER_FLAG}.dll"
            "${CMAKE_SHARED_LIBRARY_PREFIX}${LINKER_FLAG}.dll"
            "lib${LINKER_FLAG}.dll"
            "${LINKER_FLAG}.dll.a"
            "${CMAKE_SHARED_LIBRARY_PREFIX}${LINKER_FLAG}.dll.a"
            "lib${LINKER_FLAG}.dll.a")
    endif()
    unset(POSSIBLE_LIBRARY_FILENAME)
    if(IS_STATIC)
        set(TARGET_POSTFIX "STATIC")
    else()
        set(TARGET_POSTFIX "SHARED")
    endif()
    foreach(POSSIBLE_LIBRARY_FILENAME ${POSSIBLE_LIBRARY_FILENAMES})
        find_library(CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX} ${POSSIBLE_LIBRARY_FILENAME})
        if(CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX})
            message(
                DEBUG
                "${LINKER_FLAG}: Finding ${POSSIBLE_LIBRARY_FILENAME} SUCCESS = ${CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX}}!"
            )
            set(${OUTPUT_VARIABLE}
                ${CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX}}
                PARENT_SCOPE)
            unset(CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX} CACHE)
            return()
        else()
            message(DEBUG "${LINKER_FLAG}: Finding ${POSSIBLE_LIBRARY_FILENAME} FAILED!")
        endif()
    endforeach()
    set(${OUTPUT_VARIABLE}
        ${LINKER_FLAG}-NOTFOUND
        PARENT_SCOPE)
endfunction()

#[=======================================================================[
ceu_cm_get_linker_flags_from_pkg_config -- Get linker flags from pkgconfig (*.pc) files.

Synopsis: `ceu_cm_get_linker_flags_from_pkg_config(OUTPUT_VARIABLE PKGCONFIG_NAME IS_STATIC)`

Params:
    - `OUTPUT_VARIABLE`: Name of the output variable.
    - `PKGCONFIG_NAME`: pkgconfig entry. i.e., name of the `*.pc` file.
    - `IS_STATIC`: Whether to find static or dynamic libraries.

Sample:
    - `ceu_cm_get_linker_flags_from_pkg_config(OV xrender TRUE)` -> `Xrender;X11;pthread;xcb;Xau;Xdmcp`
    - `ceu_cm_get_linker_flags_from_pkg_config(OV xrender FALSE)` -> `Xrender;X11`

Sets:
    - `OUTPUT_VARIABLE`: Parent scope level.
    - `CEU_CM_PKGCONF_LIB_${PKGCONFIG_NAME}*`: Cache level. For debug purposes.
#]=======================================================================]
function(ceu_cm_get_linker_flags_from_pkg_config OUTPUT_VARIABLE PKGCONFIG_NAME IS_STATIC)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(CEU_CM_PKGCONF_LIB_${PKGCONFIG_NAME} ${PKGCONFIG_NAME} QUIET)
        if(CEU_CM_PKGCONF_LIB_${PKGCONFIG_NAME}_FOUND)
            if(NOT ${IS_STATIC} AND CEU_CM_PKGCONF_LIB_${PKGCONFIG_NAME}_LIBRARIES)
                set(${OUTPUT_VARIABLE}
                    ${CEU_CM_PKGCONF_LIB_${PKGCONFIG_NAME}_LIBRARIES}
                    PARENT_SCOPE)
                return()
            elseif(${IS_STATIC} AND CEU_CM_PKGCONF_LIB_${PKGCONFIG_NAME}_STATIC_LIBRARIES)
                set(${OUTPUT_VARIABLE}
                    ${CEU_CM_PKGCONF_LIB_${PKGCONFIG_NAME}_STATIC_LIBRARIES}
                    PARENT_SCOPE)
                return()
            endif()
        endif()
    endif()
    set(${OUTPUT_VARIABLE}
        ${OUTPUT_VARIABLE}-NOTFOUND
        PARENT_SCOPE)
endfunction()

#[=======================================================================[
ceu_cm_get_library_abspath_from_pkg_config -- Get library absolute paths from pkgconfig (*.pc) files.

Synopsis: `ceu_cm_get_library_abspath_from_pkg_config(OUTPUT_VARIABLE PKGCONFIG_NAME IS_STATIC)`

Params:
    - `OUTPUT_VARIABLE`: Name of the output variable.
    - `PKGCONFIG_NAME`: pkgconfig entry. i.e., name of the `*.pc` file.
    - `IS_STATIC`: Whether to find static or dynamic libraries.

Sample:
    - `ceu_cm_get_library_abspath_from_pkg_config(OV xrender TRUE)` -> `/usr/lib/x86_64-linux-gnu/libXrender.a;/usr/lib/x86_64-linux-gnu/libX11.a;/usr/lib/x86_64-linux-gnu/libpthread.a;/usr/lib/x86_64-linux-gnu/libxcb.a;/usr/lib/x86_64-linux-gnu/libXau.a;/usr/lib/x86_64-linux-gnu/libXdmcp.a`
    - `ceu_cm_get_library_abspath_from_pkg_config(OV xrender FALSE)` -> `/usr/lib/x86_64-linux-gnu/libXrender.so;/usr/lib/x86_64-linux-gnu/libX11.so`

Sets:
    - `OUTPUT_VARIABLE`: Parent scope level.
#]=======================================================================]
function(ceu_cm_get_library_abspath_from_pkg_config OUTPUT_VARIABLE PKGCONFIG_NAME IS_STATIC)
    ceu_cm_get_linker_flags_from_pkg_config(LINKER_FLAGS "${PKGCONFIG_NAME}" "${IS_STATIC}")
    if(NOT LINKER_FLAGS)
        set(${OUTPUT_VARIABLE}
            ${OUTPUT_VARIABLE}-NOTFOUND
            PARENT_SCOPE)
        return()
    endif()

    foreach(LINKER_FLAG ${LINKER_FLAGS})
        ceu_cm_get_abspath_from_linker_flag(${LINKER_FLAG}_LIBRARY_ABSPATH ${LINKER_FLAG} ${IS_STATIC})
        if(${LINKER_FLAG}_LIBRARY_ABSPATH)
            list(APPEND THIS_LIBRARY_ABSPATHS ${${LINKER_FLAG}_LIBRARY_ABSPATH})
        else()
            list(APPEND THIS_LIBRARY_ABSPATHS ${LINKER_FLAG}-NOTFOUND)
        endif()
    endforeach()
    set(${OUTPUT_VARIABLE}
        ${THIS_LIBRARY_ABSPATHS}
        PARENT_SCOPE)
endfunction()

#[=======================================================================[
ceu_cm_enhanced_find_library -- Get library absolute paths from pkgconfig (*.pc) files (piortized, recommended) or linker flags.

Synopsis: `ceu_cm_enhanced_find_library([STATIC] [NO_CACHE] [OUTPUT_VARIABLE OUTPUT_VARIABLE] [LINKER_FLAG LINKER_FLAG] [PKGCONFIG_NAME PKGCONFIG_NAME])`

Params:
    - `OUTPUT_VARIABLE`: Name of the output variable.
    - `PKGCONFIG_NAME`: pkgconfig entry. i.e., name of the `*.pc` file.
    - `STATIC`: If this option was set, find static libraries; otherwise find dynamic libraries.
    - `LINKER_FLAG`: Linker flag without prefixes like `-l`.
    - `NO_CACHE`:
        - If this option was set, set result to parent scope; otherwise to global cache.
        - If this option was set, do not rely on global cache. i.e., if this package had already been searched, ignore and overwrite previous results.

Sets:
    - `OUTPUT_VARIABLE`: Parent scope level.
#]=======================================================================]
function(ceu_cm_enhanced_find_library)
    set(options STATIC NO_CACHE)
    set(oneValueArgs OUTPUT_VARIABLE LINKER_FLAG PKGCONFIG_NAME)
    set(multiValueArgs "")
    cmake_parse_arguments(CEU_CM_ENHANCED_FIND_LIBRARY "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(CEU_CM_ENHANCED_FIND_LIBRARY_NO_CACHE AND DEFINED CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE)
        return()
    endif()

    if(CEU_CM_ENHANCED_FIND_LIBRARY_STATIC)
        set(CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_TYPE "static")
    else()
        set(CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_TYPE "shared")
    endif()

    if(NOT DEFINED CEU_CM_ENHANCED_FIND_LIBRARY_LINKER_FLAG AND NOT DEFINED CEU_CM_ENHANCED_FIND_LIBRARY_PKGCONFIG_NAME)
        message(FATAL_ERROR "You need to define LINKER_FLAG or LIBRARY_PKGCONFIG_NAME")
    endif()

    if(DEFINED CEU_CM_ENHANCED_FIND_LIBRARY_PKGCONFIG_NAME
       AND CEU_CM_ENHANCED_FIND_LIBRARY_PKGCONFIG_NAME
       AND PKG_CONFIG_FOUND)
        ceu_cm_get_library_abspath_from_pkg_config(
            "${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS"
            "${CEU_CM_ENHANCED_FIND_LIBRARY_PKGCONFIG_NAME}" "${CEU_CM_ENHANCED_FIND_LIBRARY_STATIC}")
    endif()
    if(${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS)
        if(CEU_CM_ENHANCED_FIND_LIBRARY_NO_CACHE)
            set(${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}
                ${${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS}
                PARENT_SCOPE)
        else()
            set(${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}
                ${${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS}
                CACHE INTERNAL
                      "Result of calling ceu_cm_enhanced_find_library(${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE})")
        endif()
        message(
            "CEU_CM: Find ${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_TYPE} ${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE} ${${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}}"
        )
        return()
    endif()

    if(DEFINED CEU_CM_ENHANCED_FIND_LIBRARY_LINKER_FLAG AND CEU_CM_ENHANCED_FIND_LIBRARY_LINKER_FLAG)
        ceu_cm_get_abspath_from_linker_flag(
            "${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS"
            "${CEU_CM_ENHANCED_FIND_LIBRARY_LINKER_FLAG}" "${CEU_CM_ENHANCED_FIND_LIBRARY_STATIC}")
    endif()
    if(NOT ${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS)
        set(${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS
            ${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}-NOTFOUND)
    endif()
    if(CEU_CM_ENHANCED_FIND_LIBRARY_NO_CACHE)
        set(${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}
            ${${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS}
            PARENT_SCOPE)
    else()
        set(${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}
            ${${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS}
            CACHE INTERNAL
                  "Result of calling ceu_cm_enhanced_find_library(${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE})")
    endif()
    message(
        "CEU_CM: Find ${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_TYPE} ${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE} ${${CEU_CM_ENHANCED_FIND_LIBRARY_OUTPUT_VARIABLE}}"
    )
    return()
endfunction()
