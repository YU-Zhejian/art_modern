#pragma once
#include <exception>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>

namespace labw {
namespace art_modern {

    struct FastaRecord {
        const std::string id;
        const std::string sequence;

        FastaRecord(const FastaRecord&) = delete;
        FastaRecord(FastaRecord&&) = default; // Allow move constructor
        FastaRecord& operator=(const FastaRecord&) = delete;
        FastaRecord& operator=(FastaRecord&&) = delete;
    };

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

} // namespace art_modern
} // namespace labw