#pragma once

#include "FileStreamInterface.hh"

class DumbFileStream : public FileStreamInterface {
public:
    explicit DumbFileStream(const std::string& filename);
    explicit DumbFileStream();

    void write(const std::string& str) override;

    void close() override;

    ~DumbFileStream() override;

private:
};
