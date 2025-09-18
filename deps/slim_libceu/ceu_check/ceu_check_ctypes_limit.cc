#include "ceu_check/ceu_check_ctypes_limit.hh"

#include <cinttypes>
#include <climits>
#include <cstdint>
#include <sstream>
#include <string>

std::string ceu_check_get_ctypes_limit_info()
{
    std::ostringstream oss;
    oss << "Compile-time C Types max, min, etc. limits:" << std::endl;
    oss << "\tchar           (" << sizeof(char) << " size):      "
        << std::showpos << static_cast<int>(CHAR_MIN) << " -> " << static_cast<int>(CHAR_MAX) << std::noshowpos << std::endl;
    oss << "\tschar          (" << sizeof(signed char) << " size):      "
        << std::showpos << static_cast<int>(SCHAR_MIN) << " -> " << static_cast<int>(SCHAR_MAX) << std::noshowpos << std::endl;
    oss << "\tuchar          (" << sizeof(unsigned char) << " size):      "
        << std::showpos << 0 << " -> " << static_cast<unsigned int>(UCHAR_MAX) << std::noshowpos << std::endl;
    oss << "\tsize_t         (" << sizeof(size_t) << " size):      "
        << std::showpos << 0 << " -> " << SIZE_MAX << std::noshowpos << std::endl;
    oss << "\tptrdiff_t      (" << sizeof(ptrdiff_t) << " size):      "
        << std::showpos << PTRDIFF_MIN << " -> " << PTRDIFF_MAX << std::noshowpos << std::endl;
    oss << "\tshort          (" << sizeof(short) << " size):      "
        << std::showpos << SHRT_MIN << " -> " << SHRT_MAX << std::noshowpos << std::endl;
    oss << "\tushort         (" << sizeof(unsigned short) << " size):      "
        << std::showpos << 0 << " -> " << USHRT_MAX << std::noshowpos << std::endl;
    oss << "\tint            (" << sizeof(int) << " size):      "
        << std::showpos << INT_MIN << " -> " << INT_MAX << std::noshowpos << std::endl;
    oss << "\tuint           (" << sizeof(unsigned int) << " size):      "
        << std::showpos << 0 << " -> " << UINT_MAX << std::noshowpos << std::endl;
    oss << "\tlong           (" << sizeof(long) << " size):      "
        << std::showpos << LONG_MIN << " -> " << LONG_MAX << std::noshowpos << std::endl;
    oss << "\tulong          (" << sizeof(unsigned long) << " size):      "
        << std::showpos << 0 << " -> " << ULONG_MAX << std::noshowpos << std::endl;
    oss << "\tllong          (" << sizeof(long long) << " size):      "
        << std::showpos << LLONG_MIN << " -> " << LLONG_MAX << std::noshowpos << std::endl;
    oss << "\tullong         (" << sizeof(unsigned long long) << " size):      "
        << std::showpos << 0 << " -> " << ULLONG_MAX << std::noshowpos << std::endl;
    oss << "\tbool           (" << sizeof(bool) << " size):      "
        << std::showpos << static_cast<int>(false) << " -> " << static_cast<int>(true) << std::noshowpos << std::endl;
    return oss.str();
}
