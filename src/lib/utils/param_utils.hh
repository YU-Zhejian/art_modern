#include "mpi_utils.hh"
#include <boost/log/trivial.hpp>
#include <boost/program_options/variables_map.hpp>

namespace labw::art_modern {

template <typename T> T get_param(const boost::program_options::variables_map& vm, const std::string& name)
{
    try {
        return vm[name].as<T>();
    } catch (const std::exception& exp) {
        BOOST_LOG_TRIVIAL(fatal) << exp.what();
        abort_mpi();
    }
}

}
