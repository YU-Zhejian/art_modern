unset(USE_STL_RANDOM)
unset(USE_BOOST_RANDOM)
unset(USE_ONEMKL_RANDOM)
unset(USE_SYSTEM_PCG_RANDOM)
unset(WITH_ONEMKL)
unset(WITH_SYSTEM_PCG)

if(NOT DEFINED USE_RANDOM_GENERATOR)
    set(USE_RANDOM_GENERATOR STL)
endif()

if("${USE_RANDOM_GENERATOR}" STREQUAL "STL")
    set(USE_STL_RANDOM ON)
elseif("${USE_RANDOM_GENERATOR}" STREQUAL "GSL")
    message(FATAL "Support for GSL random generator has been removed.")
elseif("${USE_RANDOM_GENERATOR}" STREQUAL "BOOST")
    set(USE_BOOST_RANDOM ON)
elseif("${USE_RANDOM_GENERATOR}" STREQUAL "PCG")
    set(USE_PCG_RANDOM ON)
elseif("${USE_RANDOM_GENERATOR}" STREQUAL "SYSTEM_PCG")
    set(USE_SYSTEM_PCG_RANDOM ON)
elseif("${USE_RANDOM_GENERATOR}" STREQUAL "ONEMKL")
    if(DEFINED FIND_RANDOM_MKL_THROUGH_PKGCONF)
        if(BUILD_SHARED_LIBS)
            ceu_cm_enhanced_find_library(OUTPUT_VARIABLE mkl LINKER_FLAG "${FIND_RANDOM_MKL_THROUGH_PKGCONF}"
                                         PKGCONF_NAME "${FIND_RANDOM_MKL_THROUGH_PKGCONF}")
        else()
            ceu_cm_enhanced_find_library(
                OUTPUT_VARIABLE
                mkl
                STATIC
                LINKER_FLAG
                "${FIND_RANDOM_MKL_THROUGH_PKGCONF}"
                PKGCONF_NAME
                "${FIND_RANDOM_MKL_THROUGH_PKGCONF}")
        endif()
        if(NOT TARGET CEU_CM_EFL::mkl)
            message(FATAL_ERROR "mkl (${FIND_RANDOM_MKL_THROUGH_PKGCONF}) not found!")
        endif()
        set(USE_ONEMKL_RANDOM ON)
        set(WITH_ONEMKL ON)
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} "CEU_CM_EFL::mkl")
    else()
        find_package(MKL REQUIRED)
        set(USE_ONEMKL_RANDOM ON)
        set(WITH_ONEMKL ON)
        include_directories(${MKL_INCLUDE} ${SYCL_INCLUDE_DIR})
        set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} MKL::MKL)
    endif()
else()
    message(
        FATAL_ERROR
            "Variable USE_RANDOM_GENERATOR should be set to one of STL, GSL, BOOST or ONEMKL. Current: ${USE_RANDOM_GENERATOR}"
    )
endif()
