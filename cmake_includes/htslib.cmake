if(DEFINED USE_HTSLIB)
    if(BUILD_SHARED_LIBS)
        # FIXME
        ceu_cm_enhanced_find_library(OUTPUT_VARIABLE CMAKE_CM_LIBHTS_LIBRARY LINKER_FLAG ${USE_HTSLIB} PKGCONF_NAME
                                     htslib)
    else()
        ceu_cm_enhanced_find_library(
            OUTPUT_VARIABLE
            CMAKE_CM_LIBHTS_LIBRARY
            STATIC
            LINKER_FLAG
            ${USE_HTSLIB}
            PKGCONF_NAME
            htslib)
    endif()

    if(NOT CMAKE_CM_LIBHTS_LIBRARY)
        message(FATAL_ERROR "htslib (${USE_HTSLIB}) not found!")
    endif()

    set(USED_HTSLIB_NAME "${CMAKE_CM_LIBHTS_LIBRARY}")
else()
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../deps/labw_slim_htslib")
    include_directories("${CMAKE_CURRENT_LIST_DIR}/../deps/labw_slim_htslib")
    set(USED_HTSLIB_NAME labw_slim_htslib)
endif()
set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} "${USED_HTSLIB_NAME}")
