#pragma once

#include <cstddef>
#include <exception>
#include <istream>
#include <mutex>
#include <string>
#include <utility>

namespace labw::art_modern {

using FastaRecord = std::pair<std::string, std::string>;

struct EOFException : std::exception { };
struct MalformedFastaException : std::exception {
public:
    const char* what() const noexcept override;
};

class FastaIterator {
public:
    explicit FastaIterator(std::istream& istream);

    FastaRecord next();

    FastaIterator(const FastaIterator&) = delete;
    FastaIterator(FastaIterator&&) = delete;
    FastaIterator& operator=(const FastaIterator&) = delete;
    FastaIterator& operator=(FastaIterator&&) = delete;

private:
    std::istream& _istream;
    std::size_t _lineno = 0;
    std::mutex mutex_;
};

} // namespace labw::art_modern