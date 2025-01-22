#include "libam/out/BaseFileReadOutput.hh"

#include "libam/utils/fs_utils.hh"

#include <boost/log/trivial.hpp>

#include <string>

namespace labw::art_modern {

void BaseFileReadOutput::close()
{
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' closed.";
    is_closed_ = true;
}
BaseFileReadOutput::BaseFileReadOutput(const std::string& filename)
    : filename(filename)
{
    prepare_writer(filename);
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";
}
} // namespace labw::art_modern
