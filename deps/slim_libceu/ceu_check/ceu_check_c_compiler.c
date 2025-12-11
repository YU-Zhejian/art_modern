/**
 *@file ceu_check_c_compiler.c
 * Converted from C++ implementation to C using TONGYI LINGMA.
 * MSVC secure C standard version
 */

#include "libceu_stddef.h"

#include "ceu_check/ceu_check_c_compiler.h"

#include "ceu_check/ceu_check_c_cxx_compiler_macro.h"

#include "ceu_utils/c_snprintf.h"
#include "ceu_utils/ceu_compiler_names.h"
#include "ceu_utils/ceu_constants.h" /* NOLINT: for CEU_UNDEFINED */

#include <stdlib.h>

char* ceu_check_get_c_compiler_info()
{
    char* buffer = calloc(40960, sizeof(char));
    if (buffer == NULL) {
        return NULL;
    }
    char* ptr = buffer;
    int remaining = 40960 * sizeof(char);
    int written = 0;

#if defined(CEU_REPRODUCIBLE_BUILDS) && (CEU_REPRODUCIBLE_BUILDS == 1)
    written = CEU_SNPRINTF(ptr, remaining, "C Compiled at: N/A due to reproducible build\n");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#else
    const char* date_str =
#if defined(__DATE__)
        __DATE__;
#else
        CEU_UNDEFINED;
#endif

    const char* time_str =
#if defined(__TIME__)
        __TIME__;
#else
        CEU_UNDEFINED;
#endif

    written = CEU_SNPRINTF(ptr, remaining, "C Compiled at: %s, %s\n", date_str, time_str);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

    written = CEU_SNPRINTF(ptr, remaining, "C Compiler Identification: %s\n", CEU_COMPILER_NAME);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#if defined(CEU_COMPILER_IS_TI)
    long ti_version =
#ifdef __TI_COMPILER_VERSION__
        __TI_COMPILER_VERSION__
#elif defined(__TI_COMPILER_VERSION)
        __TI_COMPILER_VERSION
#else
        0L
#endif
        ;
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %ld.%ld.%ld\n", CEU_COMPILER_NAME_TI,
        ti_version / 1000000L, (ti_version / 1000000L) % 1000L, ti_version % 1000L);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_INTEL_CLANG)
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %d.%d.%d\n",
        CEU_COMPILER_NAME_INTEL_ONEAPI_DPCPP, __INTEL_CLANG_COMPILER / 10000, __INTEL_CLANG_COMPILER % 10000 / 100,
        __INTEL_CLANG_COMPILER % 100);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_ARM_COMPILER_LINUX)
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %d.%d (%s) build %d\n",
        CEU_COMPILER_NAME_ARM_COMPILER_LINUX, __armclang_major__, __armclang_minor__, __armclang_version__,
        __ARM_LINUX_COMPILER_BUILD__);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_ARM_COMPILER_EMBEDDED)
#ifdef __ARMCC_VERSION
    // P.VV.BBBB
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %d.%d.%d\n",
        CEU_COMPILER_NAME_ARM_COMPILER_EMBEDDED, __ARMCC_VERSION / 1000000, (__ARMCC_VERSION / 10000) % 100,
        __ARMCC_VERSION % 10000);
#else
    written = CEU_SNPRINTF(
        ptr, remaining, "\t%s compatible version number: %s\n", CEU_COMPILER_NAME_ARM_COMPILER_EMBEDDED, "UNDEFINED");
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_AOCC)
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: ", CEU_COMPILER_NAME_AOCC);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }

#ifdef __aocc_major__
    written = CEU_SNPRINTF(ptr, remaining, "%d.%d.%d\n", __aocc_major__, __aocc_minor__, __aocc_patchlevel__);
#else
    written = CEU_SNPRINTF(ptr, remaining, "%s\n", CEU_UNDEFINED);
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

/** FIXME: See: <https://github.com/cpredef/predef/blob/master/Compilers.md#intel-cc> for pre-2021 and post-2021 version
 * intepretation differences. See also:
 * <https://github.com/intel/llvm/blob/6dd1bc3465612188fda216a208341869df5d7d8b/openmp/runtime/src/kmp_version.cpp#L32>
 */
#if defined(CEU_COMPILER_IS_ICC)
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %d", CEU_COMPILER_NAME_ICC, __ICC);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }

#ifdef __INTEL_COMPILER_UPDATE
    written = CEU_SNPRINTF(ptr, remaining, ".%d", __INTEL_COMPILER_UPDATE);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif
    written = CEU_SNPRINTF(ptr, remaining, "\n");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_MSVC)
    char* msc_build_ver = (char*)calloc(32, sizeof(char));
    if (msc_build_ver == NULL) {
        free(buffer);
        return NULL;
    }
#ifdef _MSC_BUILD
    CEU_SNPRINTF(msc_build_ver, 32, "%d", _MSC_BUILD);
#else
    CEU_SNPRINTF(msc_build_ver, 32, "%s", CEU_UNDEFINED);
#endif

    char* msc_internal_ver = (char*)calloc(32, sizeof(char));
    if (msc_internal_ver == NULL) {
        free(msc_build_ver);
        free(buffer);
        return NULL;
    }
#ifdef _MSC_FULL_VER
    CEU_SNPRINTF(msc_internal_ver, 32, "%d", _MSC_FULL_VER % 100000);
#else
    CEU_SNPRINTF(msc_internal_ver, 32, "%s", CEU_UNDEFINED);
#endif

    char* msc_major_ver = (char*)calloc(32, sizeof(char));
    if (msc_major_ver == NULL) {
        free(msc_internal_ver);
        free(msc_build_ver);
        free(buffer);
        return NULL;
    }
    char* msc_minor_ver = (char*)calloc(32, sizeof(char));
    if (msc_minor_ver == NULL) {
        free(msc_major_ver);
        free(msc_internal_ver);
        free(msc_build_ver);
        free(buffer);
        return NULL;
    }
#ifdef _MSC_VER
    CEU_SNPRINTF(msc_major_ver, 32, "%d", _MSC_VER / 100);
    CEU_SNPRINTF(msc_minor_ver, 32, "%d", _MSC_VER % 100);
#else
    CEU_SNPRINTF(msc_major_ver, 32, "%s", CEU_UNDEFINED);
    CEU_SNPRINTF(msc_minor_ver, 32, "%s", CEU_UNDEFINED);
#endif
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %s.%s.%s.%s, with Visual Studio ver. %s\n",
        CEU_COMPILER_NAME_MSVC, msc_major_ver, msc_minor_ver, msc_internal_ver, msc_build_ver, CEU_VISUAL_STUDIO_VER);
    free(msc_build_ver);
    free(msc_internal_ver);
    free(msc_major_ver);
    free(msc_minor_ver);

    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_NVHPC)
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %d.%d.%d\n", CEU_COMPILER_NAME_NVHPC,
        __NVCOMPILER_MAJOR__, __NVCOMPILER_MINOR__, __NVCOMPILER_PATCHLEVEL__);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_PGIC)
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %d.%d.%d\n", CEU_COMPILER_NAME_PGIC,
        __PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_TCC)
    int tcc_major = __TINYC__ / 10000;
    int tcc_minor = (__TINYC__ - tcc_major * 10000) / 100;
    int tcc_patchlevel = __TINYC__ % 100;
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %d.%d.%d\n", CEU_COMPILER_NAME_TCC,
        tcc_major, tcc_minor, tcc_patchlevel);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_EDG)
    int edg_major = __EDG_VERSION__ / 100;
    int edg_minor = (__EDG_VERSION__ - edg_major * 100);
    written = CEU_SNPRINTF(
        ptr, remaining, "\t%s compatible version number: %d.%d\n", CEU_COMPILER_NAME_EDG, edg_major, edg_minor);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_BORLAND)
    int borland_major = __BORLANDC__ / 256;
    int borland_revision = (__BORLANDC__ - 256 * borland_major) / 16 + __BORLANDC__ % 16;
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: %d.%d, with %s\n",
        CEU_COMPILER_NAME_BORLAND, borland_major, borland_revision, CEU_CPPB_VERSION);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
    written = CEU_SNPRINTF(ptr, remaining, "\t\t__CODEGEARC_VERSION__=");
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }

#ifdef __CODEGEARC_VERSION__
    written = CEU_SNPRINTF(ptr, remaining, "%X.%X.%d\n", ((__CODEGEARC_VERSION__ & 0xFF000000) >> 24),
        ((__CODEGEARC_VERSION__ & 0x00FF0000) >> 16), ((__CODEGEARC_VERSION__ & 0x0000FFFF)));
#else
    written = CEU_SNPRINTF(ptr, remaining, "=UNDEFINED\n");
#endif

    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_CLANG)
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: ", CEU_COMPILER_NAME_CLANG);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }

#ifdef __clang_major__
    written = CEU_SNPRINTF(ptr, remaining, "%d.%d.%d", __clang_major__, __clang_minor__, __clang_patchlevel__);
#else
    written = CEU_SNPRINTF(ptr, remaining, "UNKNOWN");
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
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
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(CEU_COMPILER_IS_GCC)
    written = CEU_SNPRINTF(ptr, remaining, "\t%s compatible version number: ", CEU_COMPILER_NAME_GCC);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
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
    } else {
        free(buffer);
        return NULL;
    }
#endif

#if defined(__VERSION__)
    written = CEU_SNPRINTF(ptr, remaining, "\t__VERSION__ string: %s\n", __VERSION__);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif

    return buffer;
}
