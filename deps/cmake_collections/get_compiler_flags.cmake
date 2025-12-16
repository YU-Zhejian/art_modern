if(CEU_CM_SHOULD_USE_NATIVE)
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -march=native)
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -mtune=native -mtune)
endif()

if(MSVC)
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS /Zc:__STDC__) # Define __STDC__
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS /Zc:__cplusplus) # Define __cplusplus
endif()

if(DEFINED CEU_CM_SILENT)
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -W0 -w)
else()
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -Wall)
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -Wextra)
    if(NOT MSVC)
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -pedantic -Wpedantic)
    else()
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS /permissive) # Standards conformance
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS /sdl) # Enable Additional Security Checks
    endif()
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release") # Release
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -g0)
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -O3 -O2)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo") # Release with Debug Information
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -O3 -O2)
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -g)
    ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -rdynamic)
else() # Debug, the default.
    if(NOT MSVC)
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -Og) # Add debug info
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -g3) # Add debug info
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -O0) # Stop optimization
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS -rdynamic)
    endif()

    if(MSVC)
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS /Wp64) # Visual Studio 64 bit compatibility
        ceu_cm_enhanced_check_compiler_flag(OUT_NAME CEU_CM_CC_FLAGS FLAGS /Z7) # Visual Studio
    endif()

    set(CEU_CM_IS_DEBUG
        1
        CACHE INTERNAL "") # Also set CMake variable
endif()
