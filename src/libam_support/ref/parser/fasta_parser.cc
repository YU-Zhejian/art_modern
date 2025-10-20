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

#include "fasta_parser.hh"

#include <istream>
#include <mutex>
#include <string>
#include <utility>

namespace labw::art_modern {
FastaIterator::FastaRecord FastaIterator::next()
{
    const std::scoped_lock rhs_lk(mutex_);
    std::string next_record_id;
    std::string next_record_sequence;
    std::string next_line;
    if (staged_next_record_id_.empty()) {
        // Fetch and parse the ID.
        while (true) {
            if (_istream.eof()) {
                throw EOFException();
            }
            std::getline(_istream, next_line);
            _lineno += 1;
            if (next_line.empty()) {
                continue; // Ignored
            }
            if (next_line.back() == '\r') {
                next_line.pop_back();
            }
            if (next_line[0] != '>') {
                throw MalformedFastaException("Record ID of FASTA must start with '>' at line " + std::to_string(_lineno));
            }
            // Directly extract the record ID without splitting the whole line
            next_record_id = next_line.substr(/**Exclude > **/ 1, next_line.find_first_of(" \t\f") - 1);
            if (next_record_id.empty()) {
                throw MalformedFastaException("Record ID is empty at line " + std::to_string(_lineno));
            }
            break;
        }
    } else {
        next_record_id = std::move(staged_next_record_id_);
    }
    while (true) {
        if (_istream.eof()) {
            return { std::move(next_record_id), std::move(next_record_sequence) };
        }
        std::getline(_istream, next_line);
        _lineno += 1;
        if (next_line.empty()) {
            continue; // Ignored
        }
        if (next_line.back() == '\r') {
            next_line.pop_back();
        }
        if (next_line[0] == '>') {
            // Directly extract the record ID without splitting the whole line
            staged_next_record_id_ = next_line.substr(/**Exclude > **/ 1, next_line.find_first_of(" \t\f") - 1);
            if (staged_next_record_id_.empty()) {
                throw MalformedFastaException("Record ID is empty at line " + std::to_string(_lineno));
            }
            return { std::move(next_record_id), std::move(next_record_sequence) };
        }
        next_record_sequence += next_line;
    }
}
FastaIterator::FastaIterator(std::istream& istream)
    : _istream(istream)
{
}
} // namespace labw::art_modern
