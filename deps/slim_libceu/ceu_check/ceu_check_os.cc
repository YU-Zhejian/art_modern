#include "ceu_check/ceu_check_os.hh"

#include "ceu_check/ceu_check_os_macro.h"
#include "ceu_check/ceu_constants.h"

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

#if defined(CEU_ON_WINDOWS) || defined(CEU_ON_CYGWIN)
#include <windows.h>
#endif

#include <fstream>
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
    oss << "MinGW 32-bit API: " << CEU_UNDEFINED;
#endif
#if defined(__MINGW64_VERSION_MAJOR) && defined(__MINGW64_VERSION_MINOR) && defined(__MINGW64_VERSION_BUGFIX)
    oss << ", MinGW 64-bit API: " << __MINGW64_VERSION_MAJOR << "." << __MINGW64_VERSION_MINOR << "."
        << __MINGW64_VERSION_BUGFIX;
#else
    oss << ", MinGW 64-bit API: "<< CEU_UNDEFINED;
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
    oss << CYGWIN_VERSION_DLL_MAJOR << "." << CYGWIN_VERSION_DLL_MINOR << std::endl;
#else
    oss << "\tCYGWIN API ver. " << CEU_UNDEFINED << ", with dll (" << CEU_UNDEFINED << ") ver. " << CEU_UNDEFINED << std::endl;
#endif
#endif

#ifdef CEU_ON_POSIX
    oss << "\tPOSIX.1 Version: ";
#ifdef _POSIX_VERSION
    oss << std::to_string(_POSIX_VERSION);
#else
    oss << CEU_UNDEFINED;
#endif
    oss << std::endl;
    oss << "\tPOSIX.2 Version: ";
#ifdef _POSIX2_VERSION
    oss << std::to_string(_POSIX2_VERSION);
#else
    oss << CEU_UNDEFINED;
#endif
    oss << std::endl;
    oss << "\tSingle UNIX Specification (SUS) Version: ";

#ifdef _XOPEN_VERSION
    oss << std::to_string(_XOPEN_VERSION);
#elif defined(_XOPEN_UNIX)
    oss << "unknown";
#else
    oss << CEU_UNDEFINED;
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
    struct utsname ceu_utsname{};
    uname(&ceu_utsname);
    oss << "\tPOSIX UTSINFO:" << std::endl;
    oss << "\t\tsysname=" << ceu_utsname.sysname << std::endl;
    oss << "\t\tnodename=" << ceu_utsname.nodename << std::endl;
    oss << "\t\trelease=" << ceu_utsname.release << std::endl;
    oss << "\t\tversion=" << ceu_utsname.version << std::endl;
    oss << "\t\tmachine=" << ceu_utsname.machine << std::endl;
#else
    oss << "\tPOSIX UTSINFO: " << CEU_UNDEFINED << std::endl;
#endif
#endif
#if defined(CEU_ON_WINDOWS) || defined(CEU_ON_CYGWIN)
    oss << "\tWindows OS Version Info: " << std::endl;
    OSVERSIONINFOEX osvi;
    ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));
    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
    oss << "\t\tPlatform ID: ";
    std::string platform_id {};
    if (GetVersionEx((OSVERSIONINFO*)&osvi)) {
        switch (osvi.dwPlatformId) {
        case VER_PLATFORM_WIN32_WINDOWS:
            platform_id = "Windows 95 / 98 / Me";
            break;
        case VER_PLATFORM_WIN32_NT:
            if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 0) {
                platform_id = "Windows 2000";
            } else if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 1) {
                platform_id = "Windows XP";
            } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 0) {
                platform_id = (osvi.wProductType == VER_NT_WORKSTATION ? "Windows Vista" : "Windows Server 2008");
            } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 1) {
                platform_id = (osvi.wProductType == VER_NT_WORKSTATION ? "Windows 7" : "Windows Server 2008 R2");
            } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 2) {
                platform_id = (osvi.wProductType == VER_NT_WORKSTATION ? "Windows 8" : "Windows Server 2012");
            } else if (osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 3) {
                platform_id = (osvi.wProductType == VER_NT_WORKSTATION ? "Windows 8.1" : "Windows Server 2012 R2");
            } else if (osvi.dwMajorVersion == 10 && osvi.dwMinorVersion == 0) {
                // Here is where build numbers start to matter.
                if (osvi.wProductType == VER_NT_WORKSTATION) {
                    // https://en.wikipedia.org/wiki/Windows_10_version_history
                    // https://en.wikipedia.org/wiki/Windows_11_version_history
                    if (osvi.dwBuildNumber < 10240) {
                        platform_id =  "Windows 10 before 1507 (Threshold)";
                    } else if (10240 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 10586) {
                        platform_id = "Windows 10 1507 (Threshold)";
                    } else if (10586 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 14393) {
                        platform_id = "Windows 10 1511 (Threshold 2)";
                    } else if (14393 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 15063) {
                        platform_id = "Windows 10 1607 (Redstone, Anniversary Update)";
                    } else if (15063 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 16299) {
                        platform_id = "Windows 10 1703 (Redstone 2, Creators Update)";
                    } else if (16299 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 17134) {
                        platform_id = "Windows 10 1709 (Redstone 3, Fall Creators Update)";
                    } else if (17134 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 17763) {
                        platform_id = "Windows 10 1803 (Redstone 4, April 2018 Update)";
                    } else if (17763 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 18362) {
                        platform_id = "Windows 10 1809 (Redstone 5, October 2018 Update)";
                    } else if (18362 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 18363) {
                        platform_id = "Windows 10 1903 (19H1, May 2019 Update)";
                    } else if (18363 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19041) {
                        platform_id = "Windows 10 1909 (19H2, November 2019 Update)";
                    } else if (19041 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19042) {
                        platform_id = "Windows 10 2004 (20H1, May 2020 Update)";
                    } else if (19042 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19043) {
                        platform_id = "Windows 10 20H2 (October 2020 Update)";
                    } else if (19042 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19044) {
                        platform_id = "Windows 10 21H1 (21H1, May 2021 Update)";
                    } else if (19044 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 19045) {
                        platform_id = "Windows 10 21H2 (21H2, November 2021 Update)";
                    } else if (19045 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 22000) {
                        platform_id = "Windows 10 22H2 (22H2, October 2022 Update)";
                    } else if (22000 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 22621) {
                        platform_id = "Windows 11 21H2 (Sun Valley, October 2021 Update)";
                    } else if (22621 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 22631) {
                        platform_id = "Windows 11 22H2 (Sun Valley 2, September 2022 Update)";
                    } else if (22631 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 26100) {
                        platform_id = "Windows 11 23H2 (Sun Valley 3, October 2023 Update)";
                    } else if (26100 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 26200) {
                        platform_id = "Windows 11 24H2 (Hudson Valley, October 2024 Update)";
                    } else if (26200 <= osvi.dwBuildNumber) {
                        platform_id = "Windows 11 25H2 (2025 Update)";
                    } else {
                        platform_id = "Windows 10/11 unknown";
                    }
                } else {
                    if (osvi.dwBuildNumber < 14393) {
                        platform_id = "Windows Server 2016 before 1607";
                    } else if (14393 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 17763) {
                        platform_id = "Windows Server 2016 (1607)";
                    } else if (17763 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 20348) {
                        platform_id = "Windows Server 2019 (1809)";
                    } else if (20348 <= osvi.dwBuildNumber && osvi.dwBuildNumber < 26100) {
                        platform_id = "Windows Server 2022 (20H2)";
                    } else if (26100 <= osvi.dwBuildNumber) {
                        platform_id = "Windows Server 2025 or later";
                    }
                }
            } else {
                platform_id = "Windows NT unknown";
            }
            break;
        default:
            platform_id = "Windows unknown";
            break;
        }
        oss << platform_id << std::endl;
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
    // Probe lsb-release. Should be in /etc/os-release or /usr/lib/os-release in most modern Linux distros.
#ifdef CEU_ON_GNU_LINUX
    oss << "\tLinux lsb-release info: ";
    // See which file exists
    const char* lsb_files[] = { "/etc/lsb-release", "/usr/lib/lsb-release" };
    const char* target_file = nullptr;
    for (const auto& file : lsb_files) {
        if (access(file, F_OK) != -1) {
            target_file = file;
            break;
        }
    }
    if (target_file != nullptr) {
        oss << target_file << std::endl;
        std::ifstream ifs(target_file);
        if (ifs.is_open()) {
            std::string line;
            while (std::getline(ifs, line)) {
                oss << "\t\t" << line << std::endl;
            }
            if (ifs.bad()) {
                oss << "\t\tRead failure." << std::endl;
            }
            ifs.close();
        } else {
            oss << "\t\tRead failure." << std::endl;
        }
    } else {
        oss << "not found" << std::endl;
    }
#endif
    return oss.str();
}
