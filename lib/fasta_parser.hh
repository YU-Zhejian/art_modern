#pragma once
#include <exception>
#include <sstream>
#include <string>
#include <utility>

namespace labw {
namespace art_modern {

    struct FastaRecord {
        std::string id;
        std::string sequence;
    };

    struct EOFException : std::exception {
    };
    struct MalformedFastaException : std::exception {
        MalformedFastaException(int lineno, std::string what)
            : _lineno(lineno)
            , _what(std::move(what))
        {
        }
        const char* what() const noexcept override;

    private:
        int _lineno;
        std::string _what;
    };

    class FastaIterator {
    public:
        explicit FastaIterator(std::istream& istream)
            : _istream(istream)
        {
        }

        FastaRecord next();

    private:
        std::istream& _istream;
        int _lineno = 0;
    };

}
}