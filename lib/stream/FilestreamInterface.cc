#include "stream/FileStreamInterface.hh"

#include "NotImplementedException.hh"

namespace labw {

namespace art_modern {

    void FileStreamInterface::write(const std::string& str) { throw NotImplementedException(); }

    void FileStreamInterface::close() { throw NotImplementedException(); }

    FileStreamInterface::~FileStreamInterface() = default;

}

}