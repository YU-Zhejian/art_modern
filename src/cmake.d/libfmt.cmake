if(DEFINED USE_LIBFMT)
    if(USE_LIBFMT STREQUAL "CMAKE")
        find_package(fmt REQUIRED)
        set(USED_LIBFMT_NAME "fmt::fmt (External CMake)")
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} "fmt::fmt")
    else()
        if(BUILD_SHARED_LIBS)
            ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libfmt LINKER_FLAG ${USE_LIBFMT} PKGCONF_NAME fmt)
        else()
            ceu_cm_enhanced_find_library(
                OUTPUT_VARIABLE
                libfmt
                STATIC
                LINKER_FLAG
                ${USE_LIBFMT}
                PKGCONF_NAME
                fmt)
        endif()

        if(NOT TARGET CEU_CM_EFL::libfmt)
            message(FATAL_ERROR "libfmt (${USE_LIBFMT}) not found!")
        endif()

        set(USED_LIBFMT_NAME "${USE_LIBFMT} (External)")
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} "CEU_CM_EFL::libfmt")
    endif()
else()
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../deps/slim_fmt")
    include_directories("${CMAKE_CURRENT_LIST_DIR}/../../deps/slim_fmt/include")
    set(USED_LIBFMT_NAME "slim_libfmt (Bundled)")
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} slim_libfmt)
endif()
