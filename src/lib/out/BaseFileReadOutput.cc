#include <boost/log/trivial.hpp>

#include "BaseFileReadOutput.hh"

namespace labw::art_modern {


void BaseFileReadOutput::close()
{
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' closed.";
    is_closed_ = true;
}
BaseFileReadOutput::BaseFileReadOutput(const std::string &filename): filename(filename), is_closed_(false)
{
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";
}
}