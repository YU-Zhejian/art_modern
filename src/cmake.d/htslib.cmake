if(DEFINED USE_HTSLIB)
    if(BUILD_SHARED_LIBS)
        ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libhts LINKER_FLAG ${USE_HTSLIB} PKGCONF_NAME htslib)
    else()
        ceu_cm_enhanced_find_library(
            OUTPUT_VARIABLE
            libhts
            STATIC
            LINKER_FLAG
            ${USE_HTSLIB}
            PKGCONF_NAME
            htslib)
    endif()

    if(NOT TARGET CEU_CM_EFL::libhts)
        message(FATAL_ERROR "htslib (${USE_HTSLIB}) not found!")
    endif()

    set(USED_HTSLIB_NAME "${USE_HTSLIB}")
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} "CEU_CM_EFL::libhts")
else()
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../deps/labw_slim_htslib")
    include_directories("${CMAKE_CURRENT_LIST_DIR}/../../deps/labw_slim_htslib")
    set(USED_HTSLIB_NAME labw_slim_htslib)
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} labw_slim_htslib)
endif()
