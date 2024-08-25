#include <iostream>

#include "ThreadSafeFileStream.hh"
#include <boost/log/trivial.hpp>
namespace labw {
namespace art_modern {
    void ThreadSafeFileStream::write(const std::string& str)
    {
        std::lock_guard<std::mutex> lock(mutex_);

        if (!file_.is_open()) {
            return;
        }
        file_ << str;
        file_.flush();
    }

    ThreadSafeFileStream::ThreadSafeFileStream(const std::string& filename)
        : file_(filename)
    {

        if (!file_.is_open()) {
            BOOST_LOG_TRIVIAL(fatal) << "Can not open output file: " << filename;
            exit(EXIT_FAILURE);
        }
    }

    void ThreadSafeFileStream::close()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!file_.is_open()) {
            return;
        }
        file_.close();
    }

    ThreadSafeFileStream::~ThreadSafeFileStream() { ThreadSafeFileStream::close(); }

    ThreadSafeFileStream::ThreadSafeFileStream(ThreadSafeFileStream&& instance) noexcept
    {
        std::unique_lock<std::mutex> rhs_lk(instance.mutex_);
        file_ = std::move(instance.file_);
    }
}
}