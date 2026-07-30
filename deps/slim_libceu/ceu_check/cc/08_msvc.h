#ifndef CEU_CC_MSVC_H
#define CEU_CC_MSVC_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_c_cxx_compiler_macro.h> instead!"
#endif

#if defined(_MSC_VER)
#define CEU_COMPILER_IS_MSVC
#ifndef CEU_COMPILER_NAME
#ifdef __clang__
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_MSVC_CLANG
#else
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_MSVC
#endif
#endif
#undef CEU_VISUAL_STUDIO_VER
#if _MSC_VER < 800
#define CEU_VISUAL_STUDIO_VER "unknown before 1"
#elif _MSC_VER == 800
#define CEU_VISUAL_STUDIO_VER "1"
/* 900 is Visual C++ 2.0; version 3.0 never shipped, the line went 1.0 -> 2.0 -> 4.0. */
#elif _MSC_VER == 900
#define CEU_VISUAL_STUDIO_VER "2"
#elif _MSC_VER == 1000
#define CEU_VISUAL_STUDIO_VER "4"
#elif _MSC_VER == 1100
#define CEU_VISUAL_STUDIO_VER "5"
#elif _MSC_VER == 1200
#define CEU_VISUAL_STUDIO_VER "6"
#elif _MSC_VER == 1300
#define CEU_VISUAL_STUDIO_VER "7"
#elif _MSC_VER == 1310
#define CEU_VISUAL_STUDIO_VER "7.1 (VS 2003)"
#elif _MSC_VER == 1400
#define CEU_VISUAL_STUDIO_VER "8 (VS 2005)"
#elif _MSC_VER == 1500
#define CEU_VISUAL_STUDIO_VER "9 (VS 2008)"
#elif _MSC_VER == 1600
#define CEU_VISUAL_STUDIO_VER "10 (VS 2010)"
#elif _MSC_VER == 1700
#define CEU_VISUAL_STUDIO_VER "11 (VS 2012)"
#elif _MSC_VER == 1800
#define CEU_VISUAL_STUDIO_VER "12 (VS 2013)"
#elif _MSC_VER == 1900
#define CEU_VISUAL_STUDIO_VER "14 (VS 2015)"
#elif _MSC_VER == 1910
#define CEU_VISUAL_STUDIO_VER "15 (VS 2017 Update 1 & 2)"
#elif _MSC_VER == 1911
#define CEU_VISUAL_STUDIO_VER "15 (VS 2017 Update 3 & 4)"
#elif _MSC_VER == 1912
#define CEU_VISUAL_STUDIO_VER "15 (VS 2017 Update 5)"
#elif _MSC_VER == 1913
#define CEU_VISUAL_STUDIO_VER "15 (VS 2017 Update 6)"
#elif _MSC_VER == 1914
#define CEU_VISUAL_STUDIO_VER "15 (VS 2017 Update 7)"
#elif _MSC_VER == 1915
#define CEU_VISUAL_STUDIO_VER "15 (VS 2017 Update 8)"
#elif _MSC_VER == 1916
#define CEU_VISUAL_STUDIO_VER "15 (VS 2017 Update 9)"
#elif _MSC_VER == 1920
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 RTW)"
#elif _MSC_VER == 1921
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 1)"
#elif _MSC_VER == 1922
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 2)"
#elif _MSC_VER == 1923
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 3)"
#elif _MSC_VER == 1924
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 4)"
#elif _MSC_VER == 1925
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 5)"
#elif _MSC_VER == 1926
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 6)"
#elif _MSC_VER == 1927
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 7)"
#elif _MSC_VER == 1928
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 8 & 9)"
#elif _MSC_VER == 1929
#define CEU_VISUAL_STUDIO_VER "16 (VS 2019 Update 10 & 11)"
#elif _MSC_VER == 1930
#define CEU_VISUAL_STUDIO_VER "17 (VS 2022 RTW)"
#elif _MSC_VER == 1931
#define CEU_VISUAL_STUDIO_VER "17.1 (VS 2022)"
#elif _MSC_VER == 1932
#define CEU_VISUAL_STUDIO_VER "17.2 (VS 2022 Update 2)"
#elif _MSC_VER == 1933
#define CEU_VISUAL_STUDIO_VER "17.3 (VS 2022 Update 3)"
#elif _MSC_VER == 1934
#define CEU_VISUAL_STUDIO_VER "17.4 (VS 2022 Update 4)"
#elif _MSC_VER == 1935
#define CEU_VISUAL_STUDIO_VER "17.5 (VS 2022 Update 5)"
#elif _MSC_VER == 1936
#define CEU_VISUAL_STUDIO_VER "17.6 (VS 2022 Update 6)"
#elif _MSC_VER == 1937
#define CEU_VISUAL_STUDIO_VER "17.7 (VS 2022 Update 7)"
#elif _MSC_VER == 1938
#define CEU_VISUAL_STUDIO_VER "17.8 (VS 2022 Update 8)"
#elif _MSC_VER == 1939
#define CEU_VISUAL_STUDIO_VER "17.9 (VS 2022 Update 9)"
#elif _MSC_VER == 1940
#define CEU_VISUAL_STUDIO_VER "17.10 (VS 2022 Update 10)"
#elif _MSC_VER == 1941
#define CEU_VISUAL_STUDIO_VER "17.11 (VS 2022 Update 11)"
#elif _MSC_VER == 1942
#define CEU_VISUAL_STUDIO_VER "17.12 (VS 2022 Update 12)"
#elif _MSC_VER == 1943
#define CEU_VISUAL_STUDIO_VER "17.13 (VS 2022 Update 13)"
#elif _MSC_VER == 1944
#define CEU_VISUAL_STUDIO_VER "17.14 (VS 2022 Update 14)"
/* 1941..1944 are matched by the explicit branches above, so this covers 1945..1949 only. */
#elif _MSC_VER > 1944 && _MSC_VER < 1950
#define CEU_VISUAL_STUDIO_VER "unknown later than 17.14 (VS 2022 Update 14)"
/*
 * From _MSC_VER 1950 on, MSVC versioning is decoupled from Visual Studio versioning: VS 2026 and later ship
 * several MSVC toolsets side by side (a preview, the current default, and older in-support ones), and a new
 * toolset lands every six months. One _MSC_VER therefore no longer identifies a VS product version, so these
 * branches report the MSVC Build Tools version instead -- that is what _MSC_VER actually encodes, via
 * 14.(_MSC_VER % 100). Toolset 14.50 shipped in November 2025 as a long-term release, not as an Insider build.
 */
#elif _MSC_VER == 1950
#define CEU_VISUAL_STUDIO_VER "2026 or later (MSVC Build Tools 14.50)"
#elif _MSC_VER == 1951
#define CEU_VISUAL_STUDIO_VER "2026 or later (MSVC Build Tools 14.51)"
#elif _MSC_VER == 1952
#define CEU_VISUAL_STUDIO_VER "2026 or later (MSVC Build Tools 14.52)"
#elif _MSC_VER > 1952 && _MSC_VER < 2000
#define CEU_VISUAL_STUDIO_VER "2026 or later (unknown MSVC Build Tools later than 14.52)"
#elif _MSC_VER >= 2000
#define CEU_VISUAL_STUDIO_VER "unknown (_MSC_VER major version beyond 19)"
#else
#define CEU_VISUAL_STUDIO_VER "unknown"
#endif
#endif
#endif /* CEU_CC_MSVC_H */
