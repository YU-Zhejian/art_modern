#ifndef CEU_CC_MSVC_H
#define CEU_CC_MSVC_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif

#if defined(_MSC_VER)
#define CEU_COMPILER_IS_MSVC
#ifndef CEU_COMPILER_NAME
#define CEU_COMPILER_NAME "MSVC"
#endif
#undef CEU_VISUAL_STUDIO_VER
#if _MSC_VER < 800
#define CEU_VISUAL_STUDIO_VER "unknown before 1"
#elif _MSC_VER == 800
#define CEU_VISUAL_STUDIO_VER "1"
#elif _MSC_VER == 900
#define CEU_VISUAL_STUDIO_VER "3"
#elif _MSC_VER == 1000
#define CEU_VISUAL_STUDIO_VER "4"
#elif _MSC_VER == 1100
#define CEU_VISUAL_STUDIO_VER "5"
#elif _MSC_VER == 1200
#define CEU_VISUAL_STUDIO_VER "6"
#elif _MSC_VER == 1300
#define CEU_VISUAL_STUDIO_VER "7"
#elif _MSC_VER == 1310
#define CEU_VISUAL_STUDIO_VER "7.1 (2003)"
#elif _MSC_VER == 1400
#define CEU_VISUAL_STUDIO_VER "8 (2005)"
#elif _MSC_VER == 1500
#define CEU_VISUAL_STUDIO_VER "9 (2008)"
#elif _MSC_VER == 1600
#define CEU_VISUAL_STUDIO_VER "10 (2010)"
#elif _MSC_VER == 1700
#define CEU_VISUAL_STUDIO_VER "11 (2012)"
#elif _MSC_VER == 1800
#define CEU_VISUAL_STUDIO_VER "12 (2013)"
#elif _MSC_VER == 1900
#define CEU_VISUAL_STUDIO_VER "14 (2015)"
#elif _MSC_VER == 1910
#define CEU_VISUAL_STUDIO_VER "15 (2017 Update 1 & 2)"
#elif _MSC_VER == 1911
#define CEU_VISUAL_STUDIO_VER "15 (2017 Update 3 & 4)"
#elif _MSC_VER == 1912
#define CEU_VISUAL_STUDIO_VER "15 (2017 Update 5)"
#elif _MSC_VER == 1913
#define CEU_VISUAL_STUDIO_VER "15 (2017 Update 6)"
#elif _MSC_VER == 1914
#define CEU_VISUAL_STUDIO_VER "15 (2017 Update 7)"
#elif _MSC_VER == 1915
#define CEU_VISUAL_STUDIO_VER "15 (2017 Update 8)"
#elif _MSC_VER == 1916
#define CEU_VISUAL_STUDIO_VER "15 (2017 Update 9)"
#elif _MSC_VER == 1920
#define CEU_VISUAL_STUDIO_VER "16 (2019 RTW)"
#elif _MSC_VER == 1921
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 1)"
#elif _MSC_VER == 1922
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 2)"
#elif _MSC_VER == 1923
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 3)"
#elif _MSC_VER == 1924
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 4)"
#elif _MSC_VER == 1925
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 5)"
#elif _MSC_VER == 1926
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 6)"
#elif _MSC_VER == 1927
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 7)"
#elif _MSC_VER == 1928
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 8 & 9)"
#elif _MSC_VER == 1929
#define CEU_VISUAL_STUDIO_VER "16 (2019 Update 10 & 11)"
#elif _MSC_VER == 1930
#define CEU_VISUAL_STUDIO_VER "17 (2022 RTW)"
#elif _MSC_VER == 1931
#define CEU_VISUAL_STUDIO_VER "17.1 (2022)"
#elif _MSC_VER == 1932
#define CEU_VISUAL_STUDIO_VER "17.2 (2022)"
#elif _MSC_VER == 1933
#define CEU_VISUAL_STUDIO_VER "17.3 (2022)"
#elif _MSC_VER == 1934
#define CEU_VISUAL_STUDIO_VER "17.4 (2022)"
#elif _MSC_VER == 1935
#define CEU_VISUAL_STUDIO_VER "17.5 (2022)"
#elif _MSC_VER == 1936
#define CEU_VISUAL_STUDIO_VER "17.6 (2022)"
#elif _MSC_VER == 1937
#define CEU_VISUAL_STUDIO_VER "17.7 (2022)"
#elif _MSC_VER == 1938
#define CEU_VISUAL_STUDIO_VER "17.8 (2022)"
#elif _MSC_VER == 1939
#define CEU_VISUAL_STUDIO_VER "17.9 (2022)"
#elif _MSC_VER == 1940
#define CEU_VISUAL_STUDIO_VER "17.10 (2022)"
#elif _MSC_VER > 1940
#define CEU_VISUAL_STUDIO_VER "unknown later than 17.10 (1940)"
#else
#define CEU_VISUAL_STUDIO_VER "unknown"
#endif
#endif
#endif
