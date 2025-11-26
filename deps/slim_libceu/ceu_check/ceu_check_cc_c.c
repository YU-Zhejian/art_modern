#include "ceu_check/ceu_check_cc_macro.h"
#include "libceu_stddef.h"

/**
 * Converted from C++ implementation to C using TONGYI LINGMA.
 * MSVC secure C standard version
 */

#include "ceu_check/c_snprintf.h"
#include "ceu_check/ceu_check_cc.h"
#include "libceu_stddef.h"

#include <stdlib.h>
#include <string.h>

#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x)

char* ceu_check_get_compiler_info() {
    char* buffer = calloc(40960, sizeof(char));
    if (buffer == NULL) {
        return NULL;
    }
    char* ptr = buffer;
    int remaining = 40960 * sizeof(char);
    int written;

#ifdef CEU_REPRODUCIBLE_BUILDS
    written = CEU_SNPRINTF(ptr, remaining, "Compiled at: N/A due to reproducible build\n");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#else
    const char* date_str = 
#if defined(__DATE__)
        __DATE__;
#else
        "unknown date";
#endif

    const char* time_str = 
#if defined(__TIME__)
        __TIME__;
#else
        "unknown time";
#endif

    written = CEU_SNPRINTF(ptr, remaining, "Compiled at: %s, %s\n", date_str, time_str);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

    written = CEU_SNPRINTF(ptr, remaining, "Compiler Identification:\n");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }

#if defined(CEU_COMPILER_IS_INTEL_CLANG)
    written = CEU_SNPRINTF(ptr, remaining, "\tIntel Clang compatible version number: %d.%d.%d\n",
                      __INTEL_CLANG_COMPILER / 10000,
                      __INTEL_CLANG_COMPILER % 10000 / 100,
                      __INTEL_CLANG_COMPILER % 100);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_AOCC)
    written = CEU_SNPRINTF(ptr, remaining, "\tAMD Optimizing C++ Compiler (AOCC) compatible version number: ");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
    
#ifdef __aocc_major__
    written = CEU_SNPRINTF(ptr, remaining, "%d.%d.%d\n", __aocc_major__, __aocc_minor__, __aocc_patchlevel__);
#else
    written = CEU_SNPRINTF(ptr, remaining, "unknown\n");
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_ICC)
    written = CEU_SNPRINTF(ptr, remaining, "\tIntel Compiler Classic (ICC) compatible version number: %d",
                      __ICC);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
    
#ifdef __INTEL_COMPILER_UPDATE
    written = CEU_SNPRINTF(ptr, remaining, ".%d", __INTEL_COMPILER_UPDATE);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif
    
    written = CEU_SNPRINTF(ptr, remaining, "\n");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_MSVC)
    const char* msc_build_ver = 
#ifdef _MSC_BUILD
        #define MSC_BUILD_STR_HELPER(x) #x
        #define MSC_BUILD_STR(x) MSC_BUILD_STR_HELPER(x)
        MSC_BUILD_STR(_MSC_BUILD);
#else
        "unknown";
#endif

    const char* msc_internal_ver = 
#ifdef _MSC_FULL_VER
        #define MSC_INTERNAL_STR_HELPER(x) #x
        #define MSC_INTERNAL_STR(x) MSC_INTERNAL_STR_HELPER(x)
        MSC_INTERNAL_STR(_MSC_FULL_VER % 100000);
#else
        "unknown";
#endif

    const char* msc_major_ver = 
#ifdef _MSC_VER
        #define MSC_MAJOR_STR_HELPER(x) #x
        #define MSC_MAJOR_STR(x) MSC_MAJOR_STR_HELPER(x)
        MSC_MAJOR_STR(_MSC_VER / 100);
#else
        "unknown";
#endif

    const char* msc_minor_ver = 
#ifdef _MSC_VER
        #define MSC_MINOR_STR_HELPER(x) #x
        #define MSC_MINOR_STR(x) MSC_MINOR_STR_HELPER(x)
        MSC_MINOR_STR(_MSC_VER % 100);
#else
        "unknown";
#endif

    written = CEU_SNPRINTF(ptr, remaining, "\tMSVC compatible version number: %s.%s.%s.%s, with Visual Studio ver. %s\n",
                      msc_major_ver, msc_minor_ver, msc_internal_ver, msc_build_ver, 
                      STRINGIFY(CEU_VISUAL_STUDIO_VER));
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_NVHPC)
    written = CEU_SNPRINTF(ptr, remaining, "\tNVidia High-Performance Compiler (NVHPC) compatible version number: %d.%d.%d\n",
                      __NVCOMPILER_MAJOR__, __NVCOMPILER_MINOR__, __NVCOMPILER_PATCHLEVEL__);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_PGIC)
    written = CEU_SNPRINTF(ptr, remaining, "\tPGI Compiler (PGIC) compatible version number: %d.%d.%d\n",
                      __PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_TCC)
    int tcc_major = __TINYC__ / 10000;
    int tcc_minor = (__TINYC__ - tcc_major * 10000) / 100;
    int tcc_patchlevel = __TINYC__ % 100;
    written = CEU_SNPRINTF(ptr, remaining, "\tTiny C Compiler (TCC) compatible version number: %d.%d.%d\n",
                      tcc_major, tcc_minor, tcc_patchlevel);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_EDG)
    int edg_major = __EDG_VERSION__ / 100;
    int edg_minor = (__EDG_VERSION__ - edg_major * 100);
    written = CEU_SNPRINTF(ptr, remaining, "\tEDG compatible version number: %d.%d\n",
                      edg_major, edg_minor);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_BORLAND)
    int borland_major = __BORLANDC__ / 256;
    int borland_revision = (__BORLANDC__ - 256 * borland_major) / 16 + __BORLANDC__ % 16;
    written = CEU_SNPRINTF(ptr, remaining, "\tBorland compatible version number: %d.%d, with %s\n",
                      borland_major, borland_revision, STRINGIFY(CEU_CPPB_VERSION));
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_CLANG)
    written = CEU_SNPRINTF(ptr, remaining, "\tClang compatible version number: ");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
    
#ifdef __clang_major__
    written = CEU_SNPRINTF(ptr, remaining, "%d.%d.%d",
                      __clang_major__, __clang_minor__, __clang_patchlevel__);
#else
    written = CEU_SNPRINTF(ptr, remaining, "UNKNOWN");
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
    
#ifdef __clang_version__
    written = CEU_SNPRINTF(ptr, remaining, " (%s)", __clang_version__);
#endif
    written = CEU_SNPRINTF(ptr, remaining, "\n");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_GCC)
    written = CEU_SNPRINTF(ptr, remaining, "\tGCC compatible version number: ");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        /** In case of snprintf error or truncation */
        free(buffer);
        return NULL;
    }
    
#ifdef __GNUC_PATCHLEVEL__
    written = CEU_SNPRINTF(ptr, remaining, "%d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
    written = CEU_SNPRINTF(ptr, remaining, "%d.%d\n", __GNUC__, __GNUC_MINOR__);
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

#if defined(__VERSION__)
    written = CEU_SNPRINTF(ptr, remaining, "\t__VERSION__ string: %s\n", __VERSION__);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else{
        free(buffer);
        return NULL;
    }
#endif

    return buffer;
}
