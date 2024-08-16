#pragma once

#include "FileStreamInterface.hh"
#include <fstream>
#include <mutex>
namespace labw {
namespace art_modern {
    class ThreadSafeFileStream : public FileStreamInterface {
    public:
        explicit ThreadSafeFileStream(const std::string& filename);

        void write(const std::string& str) override;

        void close() override;

        ~ThreadSafeFileStream() override;

        ThreadSafeFileStream(ThreadSafeFileStream&& instance) noexcept;

    private:
        std::ofstream file_;
        std::mutex mutex_;
    };
}
}