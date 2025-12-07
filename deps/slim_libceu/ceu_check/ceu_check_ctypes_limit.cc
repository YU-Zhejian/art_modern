#include "ceu_check/ceu_check_ctypes_limit.hh"

#include <cstddef> /* NOLINT: for size_t, ptrdiff_t */
#include <cstdlib> /* NOLINT: for NULL */
#include <iomanip>
#include <ios>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>

namespace {
template <typename T> std::string ceu_get_int_limits(const std::string& name)
{
    std::ostringstream oss;
    oss << "\t" << std::setw(12) << name << " (" << std::setw(3) << sizeof(T) << " size, " << std::setw(3)  <<sizeof(T)  * 8 << " bits): "
    << std::showpos
        << std::setw(20) << std::to_string(std::numeric_limits<T>::min()) << " -> " << std::setw(20)
        << std::to_string(std::numeric_limits<T>::max()) << std::noshowpos;
    return oss.str();
}
template <typename T> std::string ceu_get_float_limits(const std::string& name)
{
    std::ostringstream oss;
    oss << "\t" << std::setw(12) << name << " (" << std::setw(3) << sizeof(T) << " size, " << std::setw(3)  <<sizeof(T)  * 8 << " bits): "
    << std::scientific
        << std::setprecision(12) << std::setw(20) << std::numeric_limits<T>::min() << " -> " << std::setw(20)
        << std::numeric_limits<T>::max() << ", eps=" << std::setw(20) << std::numeric_limits<T>::epsilon()
        << std::defaultfloat;
    return oss.str();
}
} // namespace

std::string ceu_check_get_ctypes_limit_info()
{
    std::ostringstream oss;
    oss << "Compile-time C Types max, min, etc. limits:" << std::endl;
    oss << ceu_get_int_limits<char>("char") << std::endl;
    oss << ceu_get_int_limits<signed char>("schar") << std::endl;
    oss << ceu_get_int_limits<unsigned char>("uchar") << std::endl;
    oss << ceu_get_int_limits<size_t>("size_t") << std::endl;
    oss << ceu_get_int_limits<ptrdiff_t>("ptrdiff_t") << std::endl;
    oss << ceu_get_int_limits<short>("short") << std::endl; /* NOLINT */
    oss << ceu_get_int_limits<unsigned short>("ushort") << std::endl; /* NOLINT */
    oss << ceu_get_int_limits<int>("int") << std::endl;
    oss << ceu_get_int_limits<unsigned int>("uint") << std::endl;
    oss << ceu_get_int_limits<long>("long") << std::endl; /* NOLINT */
    oss << ceu_get_int_limits<unsigned long>("ulong") << std::endl; /* NOLINT */
    oss << ceu_get_int_limits<long long>("llong") << std::endl; /* NOLINT */
    oss << ceu_get_int_limits<unsigned long long>("ullong") << std::endl; /* NOLINT */
    oss << ceu_get_int_limits<bool>("bool") << std::endl;

    oss << ceu_get_float_limits<float>("float") << std::endl;
    oss << ceu_get_float_limits<double>("double") << std::endl;
    oss << ceu_get_float_limits<long double>("ldouble") << std::endl;
    return oss.str();
}
