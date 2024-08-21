#include "stream/FileStreamInterface.hh"

#include "ExceptionUtils.hh"
#include "NotImplementedException.hh"
#include <boost/log/trivial.hpp>

namespace labw {

namespace art_modern {

    void FileStreamInterface::write(const std::string& str)
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls FileStreamInterface::write?";
        throw_with_trace(NotImplementedException());
    }

    void FileStreamInterface::close()
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls FileStreamInterface::close?";
        throw_with_trace(NotImplementedException());
    }

    FileStreamInterface::~FileStreamInterface() = default;

}

}