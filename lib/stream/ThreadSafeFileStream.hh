#pragma once

#include "FileStreamInterface.hh"
#include <fstream>
#include <mutex>

class ThreadSafeFileStream : public FileStreamInterface {
public:
    explicit ThreadSafeFileStream(const std::string& filename);

    void write(const std::string& str) override;

    void close() override;

    ~ThreadSafeFileStream() override;
    ThreadSafeFileStream(ThreadSafeFileStream&& instance) noexcept
    {
        std::unique_lock<std::mutex> rhs_lk(instance.mutex_);
        file_ = std::move(instance.file_);
    }

private:
    std::ofstream file_;
    std::mutex mutex_;
};
