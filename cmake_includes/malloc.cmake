unset(WITH_MIMALLOC)
unset(WITH_JEMALLOC)

if(NOT DEFINED USE_MALLOC)
    set(USE_MALLOC AUTO)
endif()

if("${USE_MALLOC}" STREQUAL AUTO)
    if(BUILD_SHARED_LIBS)
        ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libjemalloc LINKER_FLAG jemalloc PKGCONF_NAME jemalloc)
    else()
        ceu_cm_enhanced_find_library(
            OUTPUT_VARIABLE
            libjemalloc
            STATIC
            LINKER_FLAG
            jemalloc
            PKGCONF_NAME
            jemalloc)
    endif()
    if(TARGET CEU_CM_EFL::libjemalloc)
        set(WITH_JEMALLOC ON)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} "CEU_CM_EFL::libjemalloc")
    else()
        find_package(mimalloc QUIET)
        if(mimalloc_FOUND)
            set(WITH_MIMALLOC ON)
            if(BUILD_SHARED_LIBS)
                set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} mimalloc)
            else()
                set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} mimalloc-static)
            endif()
        endif()
    endif()
elseif("${USE_MALLOC}" STREQUAL MIMALLOC)
    find_package(mimalloc REQUIRED)
    set(WITH_MIMALLOC ON)
    if(BUILD_SHARED_LIBS)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} mimalloc)
    else()
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} mimalloc-static)
    endif()
elseif("${USE_MALLOC}" STREQUAL JEMALLOC)
    if(BUILD_SHARED_LIBS)
        ceu_cm_enhanced_find_library(OUTPUT_VARIABLE libjemalloc LINKER_FLAG jemalloc PKGCONF_NAME jemalloc)
    else()
        ceu_cm_enhanced_find_library(
            OUTPUT_VARIABLE
            libjemalloc
            STATIC
            LINKER_FLAG
            jemalloc
            PKGCONF_NAME
            jemalloc)
    endif()
    if(TARGET CEU_CM_EFL::libjemalloc)
        set(WITH_JEMALLOC ON)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} "CEU_CM_EFL::libjemalloc")
    else()
        message(FATAL_ERROR "Required jemalloc not found!")
    endif()
elseif("${USE_MALLOC}" STREQUAL NOP)
    # Do nothing!
else()
    message(FATAL_ERROR "Unknown USE_MALLOC: ${USE_MALLOC}. Should be one of AUTO, MIMALLOC, JEMALLOC, and NOP")
endif()
