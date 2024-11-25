// FIXME: I am too tired to fix this **

#include "ceu_check/ceu_check_os.hh"
#include "libceu_stddef.h"
#include <sstream>

#ifdef CEU_ON_POSIX
// Should NOT be removed
#include <unistd.h> // NOLINT

#if defined(CEU_HAVE_INCLUDE_SYS_UTSNAME_H) && CEU_HAVE_INCLUDE_SYS_UTSNAME_H == 1
#include <sys/utsname.h>
#endif
#endif

#if defined(CEU_HAVE_INCLUDE_CYGWIN_VERSION_H) && CEU_HAVE_INCLUDE_CYGWIN_VERSION_H == 1
#include <cygwin/version.h>
#endif

#if defined(CEU_HAVE_INCLUDE__MINGW_H) && CEU_HAVE_INCLUDE__MINGW_H == 1
#include <_mingw.h>
#endif

#if defined(CEU_ON_WINDOWS)
#include <windows.h>
#endif

std::string get_compile_time_mingw_version(void)
{
#if defined(CEU_ON_CYGWIN_LIKE)
#if defined(CEU_HAVE_INCLUDE__MINGW_H) && CEU_HAVE_INCLUDE__MINGW_H == 1

    ceu_ystr_t* rets = ceu_ystr_create_from_cstr_guarantee("MinGW ver. ", 128);
#if defined(__MINGW32_MAJOR_VERSION) && defined(__MINGW32_MINOR_VERSION)
    ceu_ystr_t* mingw32_api_ver = convert_version_to_ystr2(__MINGW32_MAJOR_VERSION, __MINGW32_MINOR_VERSION);
#else
    ceu_ystr_t* mingw32_api_ver = ceu_ystr_create_from_cstr("undefined");
#endif
#if defined(__MINGW64_VERSION_MAJOR) && defined(__MINGW64_VERSION_MINOR) && defined(__MINGW64_VERSION_BUGFIX)
    ceu_ystr_t* mingw64_api_ver
        = convert_version_to_ystr3(__MINGW64_VERSION_MAJOR, __MINGW64_VERSION_MINOR, __MINGW64_VERSION_BUGFIX);
#else
    ceu_ystr_t* mingw64_api_ver = ceu_ystr_create_from_cstr("undefined");
#endif
#if defined(__MINGW64_VERSION_STATE)
    ceu_ystr_cstr_concat_inplace(mingw64_api_ver, " (");
    ceu_ystr_cstr_concat_inplace(mingw64_api_ver, __MINGW64_VERSION_STATE);
    ceu_ystr_cstr_concat_inplace(mingw64_api_ver, ")");
#endif
    ceu_ystr_cstr_concat_inplace(rets, "mingw32='");
    ceu_ystr_concat_inplace(rets, mingw32_api_ver);
    ceu_ystr_cstr_concat_inplace(rets, "'; mingw64='");
    ceu_ystr_concat_inplace(rets, mingw64_api_ver);
    ceu_ystr_cstr_concat_inplace(rets, "'");

    ceu_ystr_destroy(mingw64_api_ver);
    ceu_ystr_destroy(mingw32_api_ver);
#else
    ceu_ystr_t* rets = ceu_ystr_create_from_cstr("MinGW ver. undefined");
#endif
    return rets;
#else
    return CEU_NULL;
#endif
}

ceu_ystr_t* get_compile_time_cygwin_version(void)
{
#if defined(CEU_ON_CYGWIN_LIKE)
#if defined(CEU_HAVE_INCLUDE_CYGWIN_VERSION_H) && CEU_HAVE_INCLUDE_CYGWIN_VERSION_H == 1
    ceu_ystr_t* rets = ceu_ystr_create_from_cstr_guarantee("CYGWIN API ver. ", 128);
    ceu_ystr_t* cygwin_api_ver = convert_version_to_ystr2(CYGWIN_VERSION_API_MAJOR, CYGWIN_VERSION_API_MINOR);
    ceu_ystr_t* cygwin_dll_ver = convert_version_to_ystr2(CYGWIN_VERSION_DLL_MAJOR, CYGWIN_VERSION_DLL_MINOR);
    ceu_ystr_concat_inplace(rets, cygwin_api_ver);
    ceu_ystr_cstr_concat_inplace(rets, ", with dll (");
    ceu_ystr_cstr_concat_inplace(rets, CYGWIN_VERSION_DLL_IDENTIFIER);
    ceu_ystr_cstr_concat_inplace(rets, ") ver. ");
    ceu_ystr_concat_inplace(rets, cygwin_dll_ver);
    ceu_ystr_destroy(cygwin_dll_ver);
    ceu_ystr_destroy(cygwin_api_ver);
#else
    ceu_ystr_t* rets = ceu_ystr_create_from_cstr("CYGWIN API ver. undefined\n\tCYGWIN DLL (undefined) ver. undefined");
#endif
    return rets;
#else
    return CEU_NULL;
#endif
}

std::string get_run_time_haiku_version(void)
{
    return nullptr; // TODO
}

std::string get_run_time_windows_version(void)
{
#if defined(CEU_ON_WINDOWS) && !defined(CEU_ON_CYGWIN)
    ceu_ystr_t* readable_version_number; /// Human-readable Windows version like Microsoft Windows 10.
    ceu_ystr_t* win_ver; /// Windows version like 10.0.0.
    ceu_ystr_t* win_spver; /// Windows service pack version.

    ceu_ystr_t* rets = ceu_ystr_create_from_cstr_guarantee("Windows ver. ", 128);
    OSVERSIONINFOEX osvi;
    ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));
    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
    if (!GetVersionEx((OSVERSIONINFO*)&osvi)) {
        // GetVersionEx Failed!
        return CEU_NULL;
    }

    switch (osvi.dwPlatformId) {
    case VER_PLATFORM_WIN32_WINDOWS:
        readable_version_number = ceu_ystr_create_from_cstr("Windows 95 / 98 / Me");
        break;
    case VER_PLATFORM_WIN32_NT:
        if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 0) {
            readable_version_number = ceu_ystr_create_from_cstr("Windows 2000");
        } else if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 1) {
            readable_version_number = ceu_ystr_create_from_cstr("Windows XP");
        } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 0) {
            readable_version_number = ceu_ystr_create_from_cstr(
                osvi.wProductType == VER_NT_WORKSTATION ? "Windows Vista" : "Windows Server 2008");
        } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 1) {
            readable_version_number = ceu_ystr_create_from_cstr(
                osvi.wProductType == VER_NT_WORKSTATION ? "Windows 7" : "Windows Server 2008 R2");
        } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 2) {
            readable_version_number = ceu_ystr_create_from_cstr(
                osvi.wProductType == VER_NT_WORKSTATION ? "Windows 8" : "Windows Server 2012");
        } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 3) {
            readable_version_number = ceu_ystr_create_from_cstr(
                osvi.wProductType == VER_NT_WORKSTATION ? "Windows 8.1" : "Windows Server 2012 R2");
        } else if (osvi.dwMajorVersion == 10 && osvi.dwMinorVersion == 0) {
            readable_version_number = ceu_ystr_create_from_cstr(
                osvi.wProductType == VER_NT_WORKSTATION ? "Windows 10" : "Windows Server 2016 or later");
        } else {
            readable_version_number = ceu_ystr_create_from_cstr("Windows NT unknown");
        }
        break;
    default:
        readable_version_number = ceu_ystr_create_from_cstr("Windows unknown");
        break;
    }

    win_ver = convert_version_to_ystr3(osvi.dwMajorVersion, osvi.dwMinorVersion, osvi.dwBuildNumber);
    ceu_ystr_concat_inplace(rets, win_ver);
    ceu_ystr_cstr_concat_inplace(rets, " (");
    ceu_ystr_concat_inplace(rets, readable_version_number);
    ceu_ystr_cstr_concat_inplace(rets, ")");

    win_spver = ceu_ystr_create_from_cstr_guarantee(" Service Pack ", 128);
    if (osvi.wServicePackMajor != 0 || osvi.wServicePackMinor != 0 || osvi.szCSDVersion[0] != '\0') {
        ceu_ystr_t* win_spver_s = convert_version_to_ystr2(osvi.wServicePackMajor, osvi.wServicePackMinor);
        ceu_ystr_concat_const(win_spver, win_spver_s);
        ceu_ystr_cstr_concat_const(win_spver, " (");
        ceu_ystr_cstr_concat_const(win_spver, osvi.szCSDVersion);
        ceu_ystr_cstr_concat_const(win_spver, ")");
        ceu_ystr_destroy(win_spver_s);
    } else {
        ceu_ystr_cstr_concat_inplace(win_spver, "not installed");
    }
    ceu_ystr_concat_inplace(rets, win_spver);
    ceu_ystr_destroy(win_spver);
    ceu_ystr_destroy(win_ver);
    ceu_ystr_destroy(readable_version_number);
    return rets;
#else // Supress on CygWin.
    return nullptr;
#endif
}

ceu_ystr_t* get_run_time_posix_uts_info(void)
{
#ifdef CEU_ON_POSIX
#if defined(CEU_HAVE_INCLUDE_SYS_UTSNAME_H) && CEU_HAVE_INCLUDE_SYS_UTSNAME_H == 1
    struct utsname ceu_utsname;
    uname(&ceu_utsname);
    ceu_ystr_t* rets = ceu_ystr_create_from_cstr_guarantee("POSIX UTSINFO='", 128);
    ceu_ystr_t* utsname = ceu_ystr_cstr_join(" ", ceu_false, 5, ceu_utsname.sysname, ceu_utsname.nodename,
        ceu_utsname.release, ceu_utsname.version, ceu_utsname.machine);
    ceu_ystr_concat_inplace(rets, utsname);
    ceu_ystr_destroy(utsname);
    return rets;
#else
    return ceu_ystr_create_from_cstr("POSIX UTSINFO=undefined");
#endif
#else
    return CEU_NULL;
#endif
}

ceu_ystr_t* get_compile_time_posix_standard(void)
{

#ifdef CEU_ON_POSIX
    ceu_ystr_t* posix1_ver = ceu_ystr_create_from_cstr_guarantee("POSIX.1 Version: ", 128);
#ifdef _POSIX_VERSION
    ceu_ystr_t* posix1_ver_s = ceu_ystr_from_uint(10, _POSIX_VERSION);
#else
    ceu_ystr_t* posix1_ver_s = ceu_ystr_create_from_cstr("undefined");
#endif
    ceu_ystr_concat_inplace(posix1_ver, posix1_ver_s);
    ceu_ystr_destroy(posix1_ver_s);

    ceu_ystr_t* posix2_ver = ceu_ystr_create_from_cstr_guarantee("POSIX.2 Version: ", 128);
#ifdef _POSIX2_VERSION
    ceu_ystr_t* posix2_ver_s = ceu_ystr_from_uint(10, _POSIX2_VERSION);
#else
    ceu_ystr_t* posix2_ver_s = ceu_ystr_create_from_cstr("undefined");
#endif
    ceu_ystr_concat_inplace(posix2_ver, posix2_ver_s);
    ceu_ystr_destroy(posix2_ver_s);

    ceu_ystr_t* sus_ver = ceu_ystr_create_from_cstr_guarantee("Single UNIX Specification (SUS) Version: ", 128);
#ifdef _XOPEN_VERSION
    ceu_ystr_t* sus_ver_s = ceu_ystr_from_uint(10, _POSIX2_VERSION);
#elif defined(_XOPEN_UNIX)
    ceu_ystr_t* sus_ver_s = ceu_ystr_create_from_cstr("unknown");
#else
    ceu_ystr_t* sus_ver_s = ceu_ystr_create_from_cstr("undefined");
#endif
    ceu_ystr_concat_inplace(sus_ver, sus_ver_s);
    ceu_ystr_destroy(sus_ver_s);

    ceu_ystr_t* sep = ceu_ystr_create_from_cstr("\n\t");
    ceu_ystr_t* rets = ceu_ystr_join(sep, ceu_false, 3, posix1_ver, posix2_ver, sus_ver);
    ceu_ystr_destroy(sep);
    ceu_ystr_destroy(posix1_ver);
    ceu_ystr_destroy(posix2_ver);
    ceu_ystr_destroy(sus_ver);
    return rets;
#else
    return CEU_NULL;
#endif
}

ceu_ystr_t* ceu_check_get_compile_time_os_info(void)
{
    ceu_ystr_t* haiku_version_buff = get_run_time_haiku_version();
    ceu_ystr_t* mingw_version_buff = get_compile_time_mingw_version();
    ceu_ystr_t* cygwin_version_buff = get_compile_time_cygwin_version();
    ceu_ystr_t* posix_version_buff = get_compile_time_posix_standard();
    ceu_ystr_t* rets = ceu_ystr_create_from_cstr_guarantee("Compile-time OS info: '", 128);
    ceu_ystr_t* sep = ceu_ystr_create_from_cstr("\n\t");
    ceu_ystr_t* info;
    ceu_ystr_cstr_concat_inplace(rets, CEU_PRIMARY_OS_TYPE);
    ceu_ystr_cstr_concat_inplace(rets, "'\n\t");
    info = ceu_ystr_join(
        sep, ceu_true, 4, haiku_version_buff, mingw_version_buff, cygwin_version_buff, posix_version_buff);
    ceu_ystr_concat_inplace(rets, info);
    ceu_ystr_destroy(cygwin_version_buff);
    ceu_ystr_destroy(mingw_version_buff);
    ceu_ystr_destroy(haiku_version_buff);
    ceu_ystr_destroy(posix_version_buff);
    ceu_ystr_destroy(sep);
    ceu_ystr_destroy(info);
    return rets;
}

ceu_ystr_t* ceu_check_get_run_time_os_info(void)
{
    ceu_ystr_t* posix_uts_buff = get_run_time_posix_uts_info();
    ceu_ystr_t* windows_buff = get_run_time_windows_version();
    ceu_ystr_t* rets = ceu_ystr_create_from_cstr_guarantee("Run-time OS info: '", 128);
    ceu_ystr_t* sep = ceu_ystr_create_from_cstr("\n\t");
    ceu_ystr_t* info;
    ceu_ystr_cstr_concat_inplace(rets, CEU_PRIMARY_OS_TYPE);
    ceu_ystr_cstr_concat_inplace(rets, "'\n\t");
    info = ceu_ystr_join(sep, ceu_true, 2, posix_uts_buff, windows_buff);
    ceu_ystr_concat_inplace(rets, info);
    ceu_ystr_destroy(windows_buff);
    ceu_ystr_destroy(posix_uts_buff);
    ceu_ystr_destroy(sep);
    ceu_ystr_destroy(info);
    return rets;
}
