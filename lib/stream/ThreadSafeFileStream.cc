#include <iostream>

#include "ThreadSafeFileStream.hh"
#include "utils/mpi_utils.hh"
#include <boost/log/trivial.hpp>
namespace labw::art_modern {
void ThreadSafeFileStream::write(const std::string& str)
{
    if (!file_.is_open()) {
        return;
    }
    std::scoped_lock lock(mutex_);
    file_ << str;
}

ThreadSafeFileStream::ThreadSafeFileStream(const std::string& filename)
    : file_(filename)
{
    if (!file_.is_open()) {
        BOOST_LOG_TRIVIAL(fatal) << "Can not open output file: " << filename;
        abort_mpi();
    }
}

void ThreadSafeFileStream::close()
{
    if (!file_.is_open()) {
        return;
    }
    std::scoped_lock lock(mutex_);
    file_.flush();
    file_.close();
}

ThreadSafeFileStream::~ThreadSafeFileStream() { ThreadSafeFileStream::close(); }

ThreadSafeFileStream::ThreadSafeFileStream(ThreadSafeFileStream&& instance) noexcept
{
    std::scoped_lock rhs_lk(instance.mutex_);
    file_ = std::move(instance.file_);
}
}