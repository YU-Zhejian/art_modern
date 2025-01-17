#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/format.hpp> // NOLINT

#include "fasta_parser.hh"

#include <istream>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
FastaRecord FastaIterator::next()
{
    std::scoped_lock rhs_lk(mutex_);
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
            throw MalformedFastaException();
        }
        std::vector<std::string> parts;
        // TODO: Optimize this
        split(parts, nextLine.substr(1), boost::is_any_of(" \t\f"));
        if (!parts.empty()) {
            next_record_id = parts[0];
        } else {
            throw MalformedFastaException();
        }
        break;
    }
    while (true) {
        if (_istream.eof()) {
            return { std::move(next_record_id), std::move(next_record_sequence) };
        }
        const auto cur_pos = _istream.tellg();
        std::getline(_istream, nextLine);
        _lineno += 1;
        if (nextLine.empty()) {
            continue; // Ignored
        }
        if (nextLine.back() == '\r') {
            nextLine.pop_back();
        }
        if (nextLine[0] == '>') {
            _istream.seekg(cur_pos); // TODO: This seek may be errornous in streams
            return { std::move(next_record_id), std::move(next_record_sequence) };
        }
        next_record_sequence += nextLine;
    }
}
FastaIterator::FastaIterator(std::istream& istream)
    : _istream(istream)
{
}

const char* MalformedFastaException::what() const noexcept { return "FASTA parse error"; }
} // namespace labw::art_modern
