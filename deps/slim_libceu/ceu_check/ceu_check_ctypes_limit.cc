#include "ceu_check/ceu_check_ctypes_limit.hh"

#include <climits>
#include <cstdint>
#include <sstream>
#include <string>

#include <boost/format.hpp>

std::string ceu_check_get_ctypes_limit_info(void)
{
    std::ostringstream oss;
    oss << "Compile-time C Types max, min, etc. limits:" << std::endl;
    oss << "\t"
        << (boost::format("char           (%d bits, %zu size):      %+21d -> %+21d") % CHAR_BIT % sizeof(char)
               % CHAR_MIN % CHAR_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("schar          (%zu size):      %+21d -> %+21d") % sizeof(signed char) % SCHAR_MIN
               % SCHAR_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("uchar          (%zu size):      %+21d -> %+21d") % sizeof(unsigned char) % 0 % UCHAR_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("std::size_t    (%zu size):      %+21d -> %+21d") % sizeof(std::size_t) % 0 % SIZE_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("std::ptrdiff_t (%zu size):      %+21d -> %+21d") % sizeof(std::ptrdiff_t) % PTRDIFF_MIN
               % PTRDIFF_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("short          (%zu size):      %+21hd -> %+21hd") % sizeof(short) % SHRT_MIN % SHRT_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("ushort          (%zu size):      %+21hd -> %+21hd") % sizeof(unsigned short) % 0u
               % USHRT_MAX)
        << std::endl;
    oss << "\t" << (boost::format("int          (%zu size):      %+21d -> %+21d") % sizeof(int) % INT_MIN % INT_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("uint         (%zu size):      %+21u -> %+21u") % sizeof(unsigned int) % 0u % UINT_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("long         (%zu size):      %+21ld -> %+21ld") % sizeof(long) % LONG_MIN % LONG_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("ulong        (%zu size):      %+21ld -> %+21lu") % sizeof(unsigned long) % 0 % ULONG_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("llong        (%zu size):      %+21lld -> %+21lld") % sizeof(long long) % LLONG_MIN
               % LLONG_MAX)
        << std::endl;
    oss << "\t"
        << (boost::format("ullong       (%zu size):      %+21lld -> %+21llu") % sizeof(unsigned long long) % 0
               % ULLONG_MAX)
        << std::endl;
    oss << "\t" << (boost::format("bool        (%zu size):      %d -> %d") % sizeof(bool) % (int)false % (int)true)
        << std::endl;

    return oss.str();
}
