/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

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
    const std::scoped_lock rhs_lk(mutex_);
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
        split(parts, nextLine, boost::is_any_of(" \t\f"));
        if (!parts.empty()) {
            next_record_id = parts[0].substr(1);
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
