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

#include "libam_support/ref/batcher/FastaStreamBatcher.hh"

#include "libam_support/ds/SkipLoaderSettings.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/ref/parser/fasta_parser.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#include <cstddef>
#include <istream>
#include <limits>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
FastaStreamBatcher::FastaStreamBatcher(
    const std::size_t batch_size, std::istream& stream, const SkipLoaderSettings& sls)
    : batch_size_(batch_size)
    , fasta_iterator_(stream)
    , sls_(sls)
{
    // Skip initial lines
    for (std::size_t i = 0; i < sls_.skip_first(); ++i) {
        try {
            fasta_iterator_.next();
        } catch (EOFException&) {
            break;
        } catch (MalformedFastaException& e) {
            BOOST_LOG_TRIVIAL(fatal) << "Malformed FASTA file with error '" << e.what() << "'.";
            abort_mpi();
        }
    }
}
InMemoryFastaFetch FastaStreamBatcher::fetch()
{
    const std::scoped_lock lock(mutex_);
    std::vector<std::string> seq_names;
    std::vector<std::string> seqs;
    if (batch_size_ != std::numeric_limits<decltype(batch_size_)>::max()) {
        seq_names.reserve(batch_size_);
        seqs.reserve(batch_size_);
    }

    std::string fetch_s;
    std::string fetch_e;
    while (seq_names.size() < batch_size_) {
        try {
            auto [id, sequence] = fasta_iterator_.next();
            if (fetch_s.empty()) {
                fetch_s = id;
            }
            fetch_e = id;
            seq_names.emplace_back(std::move(id));
            seqs.emplace_back(std::move(sequence));
        } catch (EOFException&) {
            break;
        } catch (MalformedFastaException& e) {
            BOOST_LOG_TRIVIAL(fatal) << "Malformed FASTA file with error '" << e.what() << "'.";
            abort_mpi();
        }
        // Skip others
        for (std::size_t i = 0; i < sls_.skip_others(); ++i) {
            try {
                fasta_iterator_.next();
            } catch (EOFException&) {
                break;
            } catch (MalformedFastaException& e) {
                BOOST_LOG_TRIVIAL(fatal) << "Malformed FASTA file with error '" << e.what() << "'.";
                abort_mpi();
            }
        }
    }
    BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << fetch_s << " to " << fetch_e << " (" << seq_names.size()
                            << "） created";
    return { std::move(seq_names), std::move(seqs) };
}

} // namespace labw::art_modern
