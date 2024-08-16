#pragma once
#include <string>

namespace labw {
namespace art_modern {
    class FileStreamInterface {

    public:
        virtual void write(const std::string& str);
        virtual void close();

        virtual ~FileStreamInterface();
    };
}
}