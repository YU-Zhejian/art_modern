#ifndef CEU_CHECK_ARCH_MACRO_H
#define CEU_CHECK_ARCH_MACRO_H
#if defined(__x86_64__) || defined(__x86_64) || defined(_M_AMD64) || defined(__amd64__) || defined(__amd64)
#define CEU_ARCHITECTURE_X86_64
#elif defined(i386) || defined(__i386) || defined(__i386__) || defined(__i486__) || defined(__i586__)                  \
|| defined(__i686__) || defined(_M_IX86) || defined(__X86__) || defined(_X86_)
#define CEU_ARCHITECTURE_I386
#else
#define CEU_ARCHITECTURE_UNKNOWN
#endif

#endif /* CEU_CHECK_ARCH_MACRO_H */
