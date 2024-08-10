#include <iostream>

#include "ThreadSafeFileStream.hh"
#include <boost/log/trivial.hpp>

using namespace std;

void ThreadSafeFileStream::write(const std::string& str)
{
    std::lock_guard<std::mutex> lock(mutex_);
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

void ThreadSafeFileStream::close() { file_.close(); }

ThreadSafeFileStream::~ThreadSafeFileStream() { ThreadSafeFileStream::close(); }
