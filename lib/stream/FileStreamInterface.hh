#pragma once
#include <string>

namespace labw::art_modern {
class FileStreamInterface {

public:
    virtual void write(const std::string& str) = 0;
    virtual void close() = 0;

    virtual ~FileStreamInterface();
};
}
