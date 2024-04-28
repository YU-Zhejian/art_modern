#pragma once

#include "FileStreamInterface.hh"
#include <fstream>
#include <mutex>

class ThreadSafeFileStream : public FileStreamInterface {
public:
    explicit ThreadSafeFileStream(const std::string& filename);

    void write(const std::string& str) override;

    void close() override;

private:
    std::ofstream file_;
    std::mutex mutex_;
};
