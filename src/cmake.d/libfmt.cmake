if(DEFINED USE_LIBFMT)
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
            htslib)
    endif()

    if(NOT TARGET CEU_CM_EFL::libfmt)
        message(FATAL_ERROR "libfmt (${USE_LIBFMT}) not found!")
    endif()

    set(USED_LIBFMT_NAME "${USE_LIBFMT}")
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} "CEU_CM_EFL::libfmt")
else()
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../deps/slim_fmt")
    include_directories("${CMAKE_CURRENT_LIST_DIR}/../../deps/slim_fmt/include")
    set(USED_LIBFMT_NAME slim_libfmt)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} slim_libfmt)
endif()
