#include <boost/algorithm/string/regex.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>

#include "fasta_parser.hh"

namespace labw {
namespace art_modern {
    FastaRecord FastaIterator::next()
    {
        FastaRecord nextRecord;
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
            boost::trim(nextLine);
            if (nextLine[0] != '>') {
                throw MalformedFastaException(_lineno, "New name not started by >");
            } else {
                std::vector<std::string> parts;
                std::string firstPart;
                boost::split_regex(parts, nextLine.substr(1), boost::regex("\\s+"));
                if (!parts.empty()) {
                    firstPart = parts[0];
                } else {
                    throw MalformedFastaException(_lineno, "FASTA record contains empty ID.");
                }
                nextRecord.id = firstPart;
                break;
            }
        }
        while (true) {
            if (_istream.eof()) {
                return nextRecord;
            }
            auto curPos = _istream.tellg();
            std::getline(_istream, nextLine);
            _lineno += 1;
            if (nextLine.empty()) {
                continue; // Ignored
            }
            boost::trim(nextLine);
            if (nextLine[0] == '>') {
                _istream.seekg(curPos);
                return nextRecord;
            } else {
                nextRecord.sequence += nextLine;
            }
        }
    }

    const char* MalformedFastaException::what() const noexcept { return "FASTA parse error"; }
} // namespace art_modern
} // namespace labw