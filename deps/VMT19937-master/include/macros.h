#pragma once

#include <exception>
#include <iostream>
#include <sstream>

#define VERBOSE_MSG

    // store the message in str
#ifdef VERBOSE_MSG
#    define BUILDMSG(msg) \
        std::string str; \
        { \
            std::ostringstream os; \
            os  << "Error: " << __FILE__ << ":" \
                << __LINE__ << "; " \
                << msg; \
            str = os.str(); \
        }
#else
#    define BUILDMSG(msg) \
        std::string str; \
        { \
            std::ostringstream os; \
            os  << msg; \
            str = os.str(); \
        }
#endif

#define MYASSERT(cond, msg) \
    { \
        if (!(cond)) \
        { \
            BUILDMSG(msg); \
            throw std::invalid_argument(str); \
        } \
    }

#define THROW(msg) \
    { \
        BUILDMSG(msg); \
        throw std::logic_error(str); \
    }

#define NOT_IMPLEMENTED throw std::logic_error("NOT IMPLEMENTED")

// define FORCE_INLINE
#if defined(__GNUC__)
#   define FORCE_INLINE __attribute__((always_inline)) inline
#   define NO_INLINE __attribute__((noinline)) inline
#   define MAY_ALIAS __attribute__((__may_alias__))
#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)
#   define FORCE_INLINE __forceinline
#   define NO_INLINE __declspec(noinline)
#   define MAY_ALIAS
#else
#   define FORCE_INLINE inline
#   define NO_INLINE
#   define MAY_ALIAS
#endif

