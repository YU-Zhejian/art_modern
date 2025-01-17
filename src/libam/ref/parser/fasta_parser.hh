#pragma once

#include "libam/utils/class_macros_utils.hh"

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
    [[nodiscard]] const char* what() const noexcept override;
};

class FastaIterator {
public:
    explicit FastaIterator(std::istream& istream);

    FastaRecord next();

    DELETE_COPY(FastaIterator)
    DELETE_MOVE(FastaIterator)
    ~FastaIterator() = default;

private:
    std::istream& _istream;
    std::size_t _lineno = 0;
    std::mutex mutex_;
};

} // namespace labw::art_modern