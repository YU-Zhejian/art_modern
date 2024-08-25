#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/format.hpp>

#include "fasta_parser.hh"

namespace labw {
namespace art_modern {
    FastaRecord FastaIterator::next()
    {
        std::unique_lock<std::mutex> rhs_lk(mutex_);
        // FastaRecord nextRecord;
        std::string next_record_id;
        std::string next_record_sequence;
        std::string nextLine;
        while (true) {
            if (_istream.eof()) {
                throw EOFException();
            }
            std::getline(_istream, nextLine);
            _lineno += 1;
            if (nextLine.empty()) {
                continue; // Ignored
            }
            if (nextLine.back() == '\r') {
                nextLine.pop_back();
            }
            if (nextLine[0] != '>') {
                throw MalformedFastaException(_lineno, "New name not started by >");
            } else {
                std::vector<std::string> parts;
                std::string firstPart;
                boost::split(parts, nextLine.substr(1), boost::is_any_of(" \t\f"));
                if (!parts.empty()) {
                    firstPart = parts[0];
                } else {
                    throw MalformedFastaException(_lineno, "FASTA record contains empty ID.");
                }
                next_record_id = firstPart;
                break;
            }
        }
        while (true) {
            if (_istream.eof()) {
                return { next_record_id, next_record_sequence };
            }
            auto curPos = _istream.tellg();
            std::getline(_istream, nextLine);
            _lineno += 1;
            if (nextLine.empty()) {
                continue; // Ignored
            }
            if (nextLine.back() == '\r') {
                nextLine.pop_back();
            }
            if (nextLine[0] == '>') {
                _istream.seekg(curPos);
                return { next_record_id, next_record_sequence };
            } else {
                next_record_sequence += nextLine;
            }
        }
    }
    FastaIterator::FastaIterator(std::istream& istream)
        : _istream(istream)
    {
    }

    const char* MalformedFastaException::what() const noexcept { return "FASTA parse error"; }
} // namespace art_modern
} // namespace labw