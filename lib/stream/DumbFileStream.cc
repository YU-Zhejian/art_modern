#include "DumbFileStream.hh"

namespace labw::art_modern {
void DumbFileStream::write(const std::string&)
{
    // Does nothing!
}

DumbFileStream::DumbFileStream(const std::string&)
{
    // Does nothing!
}

DumbFileStream::DumbFileStream() = default;

void DumbFileStream::close() {
    // Does nothing!
};

DumbFileStream::~DumbFileStream() { DumbFileStream::close(); }
}
