#include <iostream>

#include "DumbFileStream.hh"

using namespace std;
namespace labw {
namespace art_modern {
    void DumbFileStream::write(const std::string& str)
    {
        // Does nothing!
    }

    DumbFileStream::DumbFileStream(const std::string& filename)
    {
        // Does nothing!
    }

    DumbFileStream::DumbFileStream() = default;

    void DumbFileStream::close() {
        // Does nothing!
    };

    DumbFileStream::~DumbFileStream() { DumbFileStream::close(); }
}
}