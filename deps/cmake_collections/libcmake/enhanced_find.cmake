# Finds packages.

find_package(PkgConfig QUIET)

#[=======================================================================[
ceu_cm_get_abspath_from_linker_flag -- Get absolute path of libraries from linker flag.

Synopsis: `ceu_cm_get_abspath_from_linker_flag(OUTPUT_VARIABLE LINKER_FLAG IS_STATIC)`

Params:
    - `OUTPUT_VARIABLE`: Name of the output variable.
    - `LINKER_FLAG`: Linker flag without prefixes like `-l`.
    - `IS_STATIC`: Whether to find static or dynamic libraries.

Sample:
    - `ceu_cm_get_abspath_from_linker_flag(OV z ON)` -> `libz.a`
    - `ceu_cm_get_abspath_from_linker_flag(OV z OFF)` -> `libz.so`

Sets:
    - `OUTPUT_VARIABLE`: Parent scope level.
    - `CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${IS_STATIC}`: Cache level. For debug purposes.

Warnings:
    - Finding libraries using linker flag cannot guarantee usability!

        For example, we have `app` -> `liba.so` -> `libb.so`.

        If `liba.so` does not contain implementations defined in `libb.so`,
        direct execution of `app` would yield bugs.
#]=======================================================================]
function(ceu_cm_get_abspath_from_linker_flag OUTPUT_VARIABLE LINKER_FLAG IS_STATIC LIBDIRS)
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
            "${LINKER_FLAG}"
            "${CMAKE_SHARED_LIBRARY_PREFIX}${LINKER_FLAG}"
            "${LINKER_FLAG}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${CMAKE_SHARED_LIBRARY_PREFIX}${LINKER_FLAG}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "lib${LINKER_FLAG}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            "${LINKER_FLAG}.so"
            "${CMAKE_SHARED_LIBRARY_PREFIX}${LINKER_FLAG}.so"
            "lib${LINKER_FLAG}.so"
            "lib${LINKER_FLAG}.dylib"
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
        find_library(CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX} ${POSSIBLE_LIBRARY_FILENAME}
                     HINTS "${LIBDIRS}")
        if(CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX})
            message(
                DEBUG
                "${LINKER_FLAG}: Finding ${POSSIBLE_LIBRARY_FILENAME} SUCCESS = ${CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX}}!"
            )
            set(${OUTPUT_VARIABLE}
                ${CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX}}
                PARENT_SCOPE)
            message(
                DEBUG
                "CEU_CM_EFL: Searching: -l${LINKER_FLAG} -> ${CEU_CM_LIB_FL_${LINKER_FLAG}_ABSPATH_${TARGET_POSTFIX}} (${TARGET_POSTFIX})"
            )
            unset(POSSIBLE_LIBRARY_FILENAMES)
            return()
        else()
            message(DEBUG "${LINKER_FLAG}: Finding ${POSSIBLE_LIBRARY_FILENAME} FAILED!")
        endif()
    endforeach()
    set(${OUTPUT_VARIABLE}
        ${LINKER_FLAG}-NOTFOUND
        PARENT_SCOPE)
    unset(POSSIBLE_LIBRARY_FILENAMES)
    message(DEBUG "CEU_CM_EFL: Searching: -l${LINKER_FLAG} -> ${LINKER_FLAG}-NOTFOUND (${TARGET_POSTFIX})")
endfunction()

#[=======================================================================[
ceu_cm_get_linker_flags_from_pkg_config -- Get linker flags from pkgconfig (*.pc) files.

Synopsis: `ceu_cm_get_linker_flags_from_pkg_config(OUTPUT_VARIABLE PKGCONF_NAME IS_STATIC)`

Params:
    - `OUTPUT_VARIABLE`: Name of the output variable.
    - `PKGCONF_NAME`: pkgconfig entry. i.e., name of the `*.pc` file.
    - `IS_STATIC`: Whether to find static or dynamic libraries.

Sample:
    - `ceu_cm_get_linker_flags_from_pkg_config(OV xrender ON)` -> `Xrender;X11;pthread;xcb;Xau;Xdmcp`
    - `ceu_cm_get_linker_flags_from_pkg_config(OV xrender OFF)` -> `Xrender;X11`

Sets:
    - `OUTPUT_VARIABLE`: Parent scope level.
    - `CEU_CM_PKGCONF_LIB_${PKGCONF_NAME}*`: Cache level. For debug purposes.

Assumes:
    - `PKG_CONFIG_FOUND` is ON.
#]=======================================================================]
function(ceu_cm_get_linker_flags_from_pkg_config OUTPUT_VARIABLE PKGCONF_NAME IS_STATIC)
    pkg_check_modules(CEU_CM_PKGCONF_LIB_${PKGCONF_NAME} ${PKGCONF_NAME} QUIET)
    if(CEU_CM_PKGCONF_LIB_${PKGCONF_NAME}_FOUND)
        if(NOT ${IS_STATIC} AND CEU_CM_PKGCONF_LIB_${PKGCONF_NAME}_LIBRARIES)
            set(${OUTPUT_VARIABLE} ${CEU_CM_PKGCONF_LIB_${PKGCONF_NAME}_LIBRARIES})
        elseif(${IS_STATIC} AND CEU_CM_PKGCONF_LIB_${PKGCONF_NAME}_STATIC_LIBRARIES)
            set(${OUTPUT_VARIABLE} ${CEU_CM_PKGCONF_LIB_${PKGCONF_NAME}_STATIC_LIBRARIES})
        else()
            set(${OUTPUT_VARIABLE} ${OUTPUT_VARIABLE}-NOTFOUND)
        endif()
    else()
        set(${OUTPUT_VARIABLE} ${OUTPUT_VARIABLE}-NOTFOUND)
    endif()
    message(DEBUG "CEU_CM_EFL: Searching pkgconf ${PKGCONF_NAME} (STATIC=${IS_STATIC}) -> ${${OUTPUT_VARIABLE}}")
    set(${OUTPUT_VARIABLE}
        ${${OUTPUT_VARIABLE}}
        PARENT_SCOPE)
    unset(PKGCONF_NAME)
    unset(OUTPUT_VARIABLE)
    unset(IS_STATIC)
    return()
endfunction()

#[=======================================================================[
ceu_cm_get_library_abspath_from_pkg_config -- Get library absolute paths from pkgconfig (*.pc) files.

Synopsis: `ceu_cm_get_library_abspath_from_pkg_config(OUTPUT_VARIABLE PKGCONF_NAME IS_STATIC)`

Params:
    - `OUTPUT_VARIABLE`: Name of the output variable.
    - `PKGCONF_NAME`: pkgconfig entry. i.e., name of the `*.pc` file.
    - `IS_STATIC`: Whether to find static or dynamic libraries.

Sample:
    - `ceu_cm_get_library_abspath_from_pkg_config(OV xrender ON)` -> `/usr/lib/x86_64-linux-gnu/libXrender.a;/usr/lib/x86_64-linux-gnu/libX11.a;/usr/lib/x86_64-linux-gnu/libpthread.a;/usr/lib/x86_64-linux-gnu/libxcb.a;/usr/lib/x86_64-linux-gnu/libXau.a;/usr/lib/x86_64-linux-gnu/libXdmcp.a`
    - `ceu_cm_get_library_abspath_from_pkg_config(OV xrender OFF)` -> `/usr/lib/x86_64-linux-gnu/libXrender.so;/usr/lib/x86_64-linux-gnu/libX11.so`

Sets:
    - `OUTPUT_VARIABLE`: Parent scope level.

Assumes:
    - `PKG_CONFIG_FOUND` is ON.
#]=======================================================================]
function(ceu_cm_get_library_abspath_from_pkg_config OUTPUT_VARIABLE PKGCONF_NAME IS_STATIC)
    ceu_cm_get_linker_flags_from_pkg_config(LINKER_FLAGS "${PKGCONF_NAME}" "${IS_STATIC}")
    if(NOT LINKER_FLAGS)
        set(${OUTPUT_VARIABLE}
            ${OUTPUT_VARIABLE}-NOTFOUND
            PARENT_SCOPE)
        unset(LINKER_FLAGS)
        return()
    endif()

    foreach(LINKER_FLAG ${LINKER_FLAGS})
        if(IS_STATIC)
            set(THIS_LIBDIRS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_STATIC_LIBRARY_DIRS})
        else()
            set(THIS_LIBDIRS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_LIBRARY_DIRS})
        endif()
        ceu_cm_get_abspath_from_linker_flag(${LINKER_FLAG}_LIBRARY_ABSPATH ${LINKER_FLAG} ${IS_STATIC} ${THIS_LIBDIRS})
        if(${LINKER_FLAG}_LIBRARY_ABSPATH)
            list(APPEND THIS_LIBRARY_ABSPATHS ${${LINKER_FLAG}_LIBRARY_ABSPATH})
        else()
            message(
                STATUS
                    "CEU_CM_EFL: Searching pkgconf ${PKGCONF_NAME} (STATIC=${IS_STATIC}): -l${LINKER_FLAG} -> NOTFOUND")
            set(${OUTPUT_VARIABLE}
                ${OUTPUT_VARIABLE}-NOTFOUND
                PARENT_SCOPE)
            unset(THIS_LIBDIRS)
            return()
        endif()
    endforeach()
    unset(THIS_LIBDIRS)
    set(${OUTPUT_VARIABLE}
        ${THIS_LIBRARY_ABSPATHS}
        PARENT_SCOPE)
    return()
endfunction()

#[=======================================================================[
ceu_cm_enhanced_find_library -- Get library absolute paths from pkgconfig (*.pc) files (piortized, recommended) or linker flags.

Synopsis: `ceu_cm_enhanced_find_library([STATIC] [OUTPUT_VARIABLE OUTPUT_VARIABLE] [LINKER_FLAG LINKER_FLAG] [PKGCONF_NAME PKGCONF_NAME])`

Params:
    - `OUTPUT_VARIABLE`: Name of the output variable.
    - `PKGCONF_NAME`: pkgconfig entry. i.e., name of the `*.pc` file.
    - `STATIC`: If this option was set, find static libraries; otherwise find dynamic libraries.
    - `LINKER_FLAG`: Linker flag without prefixes like `-l`.

Sets:
    - Imported library CEU_CM_EFL::${CEU_CM_EFL_OUTPUT_VARIABLE}.
#]=======================================================================]
function(ceu_cm_enhanced_find_library)
    set(options STATIC)
    set(oneValueArgs OUTPUT_VARIABLE LINKER_FLAG PKGCONF_NAME)
    set(multiValueArgs "")
    cmake_parse_arguments(CEU_CM_EFL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(CEU_CM_EFL_STATIC)
        set(CEU_CM_EFL_OUTPUT_TYPE "STATIC")
    else()
        set(CEU_CM_EFL_OUTPUT_TYPE "SHARED")
    endif()

    if(NOT DEFINED CEU_CM_EFL_LINKER_FLAG AND NOT DEFINED CEU_CM_EFL_PKGCONF_NAME)
        message(FATAL_ERROR "You need to define LINKER_FLAG or LIBRARY_PKGCONF_NAME")
    endif()

    if(DEFINED CEU_CM_EFL_PKGCONF_NAME
       AND CEU_CM_EFL_PKGCONF_NAME
       AND PKG_CONFIG_FOUND)
        ceu_cm_get_library_abspath_from_pkg_config("${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS"
                                                   "${CEU_CM_EFL_PKGCONF_NAME}" "${CEU_CM_EFL_STATIC}")
    endif()
    if(DEFINED ${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS
       AND ${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS)
        if(CEU_CM_EFL_STATIC)
            set(THIS_CFLAGS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_STATIC_CFLAGS_OTHER})
            set(THIS_LDFLAGS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_STATIC_LDFLAGS_OTHER})
            set(THIS_LIBDIRS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_STATIC_LIBRARY_DIRS})
            set(THIS_INCDIRS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_STATIC_INCLUDE_DIRS})
        else()
            set(THIS_CFLAGS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_CFLAGS_OTHER})
            set(THIS_LDFLAGS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_LDFLAGS_OTHER})
            set(THIS_LIBDIRS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_LIBRARY_DIRS})
            set(THIS_INCDIRS ${CEU_CM_PKGCONF_LIB_${CEU_CM_EFL_PKGCONF_NAME}_INCLUDE_DIRS})
        endif()
        message(
            STATUS
                "CEU_CM_EFL: Exporting target CEU_CM_EFL::${CEU_CM_EFL_OUTPUT_VARIABLE} (${CEU_CM_EFL_OUTPUT_TYPE}): LIBS=${${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS} CFLAGS=${THIS_CFLAGS} LDFLAGS=${THIS_LDFLAGS} LIBDIRS=${THIS_LIBDIRS} INCDIRS=${THIS_INCDIRS}"
        )
        add_library(CEU_CM_EFL::${CEU_CM_EFL_OUTPUT_VARIABLE} INTERFACE IMPORTED)
        set_target_properties(
            CEU_CM_EFL::${CEU_CM_EFL_OUTPUT_VARIABLE}
            PROPERTIES INTERFACE_LINK_LIBRARIES "${${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS}"
                       INTERFACE_COMPILE_OPTIONS "${THIS_CFLAGS}"
                       INTERFACE_LINK_OPTIONS "${THIS_LDFLAGS}"
                       INTERFACE_INCLUDE_DIRECTORIES "${THIS_INCDIRS}")
        unset(THIS_CFLAGS)
        unset(THIS_LDFLAGS)
        unset(THIS_LIBDIRS)
        unset(THIS_INCDIRS)

        return()
    endif()

    if(DEFINED CEU_CM_EFL_LINKER_FLAG AND CEU_CM_EFL_LINKER_FLAG)
        ceu_cm_get_abspath_from_linker_flag("${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS"
                                            "${CEU_CM_EFL_LINKER_FLAG}" "${CEU_CM_EFL_STATIC}" "")
    endif()
    if(NOT ${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS)
        set(${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS ${CEU_CM_EFL_OUTPUT_VARIABLE}-NOTFOUND)

        unset(CEU_CM_EFL_OUTPUT_TYPE)
        unset(CEU_CM_EFL_STATIC)
        unset(CEU_CM_EFL_OUTPUT_VARIABLE)
        unset(CEU_CM_EFL_PKGCONF_NAME)
        unset(CEU_CM_EFL_LINKER_FLAG)
        return()
    endif()
    add_library(CEU_CM_EFL::${CEU_CM_EFL_OUTPUT_VARIABLE} ${CEU_CM_EFL_OUTPUT_TYPE} IMPORTED)
    set_target_properties(CEU_CM_EFL::${CEU_CM_EFL_OUTPUT_VARIABLE}
                          PROPERTIES IMPORTED_LOCATION "${${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS}")

    message(
        STATUS
            "CEU_CM_EFL: Exporting target CEU_CM_EFL::${CEU_CM_EFL_OUTPUT_VARIABLE}(${CEU_CM_EFL_OUTPUT_TYPE}) LIBS=${${CEU_CM_EFL_OUTPUT_VARIABLE}_TMP_LIBRARY_ABSPATHS} "
    )
    unset(CEU_CM_EFL_OUTPUT_TYPE)
    unset(CEU_CM_EFL_STATIC)
    unset(CEU_CM_EFL_OUTPUT_VARIABLE)
    unset(CEU_CM_EFL_PKGCONF_NAME)
    unset(CEU_CM_EFL_LINKER_FLAG)
    return()
endfunction()
