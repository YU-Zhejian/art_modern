#pragma once

#include "FileStreamInterface.hh"
namespace labw {
namespace art_modern {
    class DumbFileStream : public FileStreamInterface {
    public:
        explicit DumbFileStream(const std::string& filename);
        explicit DumbFileStream();

        void write(const std::string& str) override;

        void close() override;

        ~DumbFileStream() override;
    };
}
}