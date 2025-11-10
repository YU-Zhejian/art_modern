#include "ceu_check/ceu_check_os.hh"
#include "libceu_stddef.h"

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

#include <ostream>
#include <sstream>
#include <string>

std::string ceu_check_get_compile_time_os_info()
{
    std::stringstream oss;
    oss << "Compile-time OS info:" << std::endl;
    oss << "\tPRIMARY OS='" << CEU_PRIMARY_OS_TYPE << "'" << std::endl;
#if defined(CEU_ON_MINGW32) || defined(CEU_ON_MINGW64)
#if defined(CEU_HAVE_INCLUDE__MINGW_H) && CEU_HAVE_INCLUDE__MINGW_H == 1
    oss << "\tMinGW ver. ";
#if defined(__MINGW32_MAJOR_VERSION) && defined(__MINGW32_MINOR_VERSION)
    oss << "MinGW 32-bit API: " << __MINGW32_MAJOR_VERSION << "." << __MINGW32_MINOR_VERSION;
#else
    oss << "MinGW 32-bit API: undefined";
#endif
#if defined(__MINGW64_VERSION_MAJOR) && defined(__MINGW64_VERSION_MINOR) && defined(__MINGW64_VERSION_BUGFIX)
    oss << ", MinGW 64-bit API: " << __MINGW64_VERSION_MAJOR << "." << __MINGW64_VERSION_MINOR << "."
        << __MINGW64_VERSION_BUGFIX;
#else
    oss << ", MinGW 64-bit API: undefined";
#endif
#if defined(__MINGW64_VERSION_STATE)
    oss << " (" << __MINGW64_VERSION_STATE << ")";
#endif
    oss << std::endl;
#endif
#endif
#if defined(CEU_ON_CYGWIN)
#if defined(CEU_HAVE_INCLUDE_CYGWIN_VERSION_H) && CEU_HAVE_INCLUDE_CYGWIN_VERSION_H == 1
    oss << "\tCYGWIN API ver. ";
    oss << CYGWIN_VERSION_API_MAJOR << "." << CYGWIN_VERSION_API_MINOR;
    oss << ", with dll (" << CYGWIN_VERSION_DLL_IDENTIFIER << ") ver. ";
    oss << CYGWIN_VERSION_DLL_MAJOR << "." << CYGWIN_VERSION_DLL_MINOR;
#else
    oss << "\tCYGWIN API ver. undefined\n\tCYGWIN DLL (undefined) ver. undefined";
#endif
#endif

#ifdef CEU_ON_POSIX
    oss << "\tPOSIX.1 Version: ";
#ifdef _POSIX_VERSION
    oss << std::to_string(_POSIX_VERSION);
#else
    oss << "undefined";
#endif
    oss << std::endl;
    oss << "\tPOSIX.2 Version: ";
#ifdef _POSIX2_VERSION
    oss << std::to_string(_POSIX2_VERSION);
#else
    oss << "undefined";
#endif
    oss << std::endl;
    oss << "\tSingle UNIX Specification (SUS) Version: ";

#ifdef _XOPEN_VERSION
    oss << std::to_string(_XOPEN_VERSION);
#elif defined(_XOPEN_UNIX)
    oss << "unknown";
#else
    oss << "undefined";
#endif
    oss << std::endl;
#endif
    return oss.str();
}

std::string ceu_check_get_run_time_os_info()
{
    std::stringstream oss;
    oss << "Run-time OS info:" << std::endl;
#ifdef CEU_ON_POSIX
#if defined(CEU_HAVE_INCLUDE_SYS_UTSNAME_H) && CEU_HAVE_INCLUDE_SYS_UTSNAME_H == 1
    struct utsname ceu_utsname;
    uname(&ceu_utsname);
    oss << "\tPOSIX UTSINFO:" << std::endl;
    oss << "\t\tsysname=" << ceu_utsname.sysname << std::endl;
    oss << "\t\tnodename=" << ceu_utsname.nodename << std::endl;
    oss << "\t\trelease=" << ceu_utsname.release << std::endl;
    oss << "\t\tversion=" << ceu_utsname.version << std::endl;
    oss << "\t\tmachine=" << ceu_utsname.machine << std::endl;
#else
    oss << "\tPOSIX UTSINFO: undefined" << std::endl;
#endif
#endif
#if defined(CEU_ON_WINDOWS) && !defined(CEU_ON_CYGWIN)

    oss << "\tWindows OS Version Info: " << std::endl;
    OSVERSIONINFOEX osvi;
    ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));
    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
    oss << "\t\tPlatform ID: ";
    if (GetVersionEx((OSVERSIONINFO*)&osvi)) {
        switch (osvi.dwPlatformId) {
        case VER_PLATFORM_WIN32_WINDOWS:
            oss << "Windows 95 / 98 / Me" << std::endl;
            break;
        case VER_PLATFORM_WIN32_NT:
            if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 0) {
                oss << "Windows 2000" << std::endl;
            } else if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 1) {
                oss << "Windows XP" << std::endl;
            } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 0) {
                oss << (osvi.wProductType == VER_NT_WORKSTATION ? "Windows Vista" : "Windows Server 2008") << std::endl;
            } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 1) {
                oss << (osvi.wProductType == VER_NT_WORKSTATION ? "Windows 7" : "Windows Server 2008 R2") << std::endl;
            } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 2) {
                oss << (osvi.wProductType == VER_NT_WORKSTATION ? "Windows 8" : "Windows Server 2012") << std::endl;
            } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 3) {
                oss << (osvi.wProductType == VER_NT_WORKSTATION ? "Windows 8.1" : "Windows Server 2012 R2")
                    << std::endl;
            } else if (osvi.dwMajorVersion == 10 && osvi.dwMinorVersion == 0) {
                // Here is where build numbers start to matter.
                if (osvi.wProductType == VER_NT_WORKSTATION) {
                    // https://en.wikipedia.org/wiki/Windows_10_version_history
                    // https://en.wikipedia.org/wiki/Windows_11_version_history
                    if (osvi.dwBuildNumber < 10240) {
                        oss << "Windows 10 before 1507 (Threshold)" << std::endl;
                    } else if (10240 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 10586) {
                        oss << "Windows 10 1507 (Threshold)" << std::endl;
                    } else if (10586 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 14393) {
                        oss << "Windows 10 1511 (Threshold 2)" << std::endl;
                    } else if (14393 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 15063) {
                        oss << "Windows 10 1607 (Redstone, Anniversary Update)" << std::endl;
                    } else if (15063 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 16299) {
                        oss << "Windows 10 1703 (Redstone 2, Creators Update)" << std::endl;
                    } else if (16299 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 17134) {
                        oss << "Windows 10 1709 (Redstone 3, Fall Creators Update)" << std::endl;
                    } else if (17134 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 17763) {
                        oss << "Windows 10 1803 (Redstone 4, April 2018 Update)" << std::endl;
                    } else if (17763 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 18362) {
                        oss << "Windows 10 1809 (Redstone 5, October 2018 Update)" << std::endl;
                    } else if (18362 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 18363) {
                        oss << "Windows 10 1903 (19H1, May 2019 Update)" << std::endl;
                    } else if (18363 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19041) {
                        oss << "Windows 10 1909 (19H2, November 2019 Update)" << std::endl;
                    } else if (19041 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19042) {
                        oss << "Windows 10 2004 (20H1, May 2020 Update)" << std::endl;
                    } else if (19042 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19043) {
                        oss << "Windows 10 20H2 (October 2020 Update)" << std::endl;
                    } else if (19042 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19044) {
                        oss << "Windows 10 21H1 (21H1, May 2021 Update)" << std::endl;
                    } else if (19044 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19045) {
                        oss << "Windows 10 21H2 (21H2, November 2021 Update)" << std::endl;
                    } else if (19045 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 22000) {
                        oss << "Windows 10 22H2 (22H2, October 2022 Update)" << std::endl;
                    } else if (22000 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 22621) {
                        oss << "Windows 11 21H2 (Sun Valley, October 2021 Update)" << std::endl;
                    } else if (22621 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 22631) {
                        oss << "Windows 11 22H2 (Sun Valley 2, September 2022 Update)" << std::endl;
                    } else if (22631 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 26100) {
                        oss << "Windows 11 23H2 (Sun Valley 3, October 2023 Update)" << std::endl;
                    } else if (26100 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 26200) {
                        oss << "Windows 11 24H2 (Hudson Valley, October 2024 Update)" << std::endl;
                    } else if (26200 <= osvi.dwBuildNumber) {
                        oss << "Windows 11 25H2 (2025 Update)" << std::endl;
                    } else {
                        oss << "Windows 10/11 unknown" << std::endl;
                    }
                } else {
                    if (osvi.dwBuildNumber < 14393){
                        oss << "Windows Server 2016 before 1607" << std::endl;
                    } else if (14393 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 17763) {
                        oss << "Windows Server 2016 (1607)" << std::endl;
                    } else if (17763 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 20348) {
                        oss << "Windows Server 2019 (1809)" << std::endl;
                    } else if (20348 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 26100) {
                        oss << "Windows Server 2022 (20H2)" << std::endl;
                    } else if (26100 <= osvi.dwBuildNumber) {
                        oss << "Windows Server 2025 or later" << std::endl;
                    }
                }
            } else {
                oss << "Windows NT unknown" << std::endl;
            }
            break;
        default:
            oss << "Windows unknown" << std::endl;
            break;
        }
        oss << "\t\tMajor Version: " << osvi.dwMajorVersion << "." << osvi.dwMinorVersion << "." << osvi.dwBuildNumber
            << std::endl;
        oss << "\t\tService Pack: ";
        if (osvi.wServicePackMajor != 0 || osvi.wServicePackMinor != 0 || osvi.szCSDVersion[0] != '\0') {
            oss << osvi.wServicePackMajor << "." << osvi.wServicePackMinor;
            oss << " (" << osvi.szCSDVersion << ")" << std::endl;
        } else {
            oss << "not installed" << std::endl;
        }
    } else {
        oss << "failed." << std::endl;
    }

#endif
    return oss.str();
}
