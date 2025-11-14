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
#include <art_modern_config.h> // NOLINT: For CEU_CM_IS_DEBUG

#include "libam_support/ref/fetch/FaidxFetch.hh"

#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <fmt/format.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>

#include <htslib/faidx.h>
#include <htslib/hts.h>

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace labw::art_modern {
namespace {

    std::vector<std::string> get_seq_names(const faidx_t* faidx)
    {
        std::vector<std::string> seq_names;
        const auto size = faidx_nseq(faidx);
        seq_names.reserve(size);

        for (int i = 0; i < size; i++) {
            const auto* const seq_name = faidx_iseq(faidx, i);
            if (seq_name == nullptr) {
                BOOST_LOG_TRIVIAL(fatal) << "Sequence name of seq " << i << " is null!";
                abort_mpi();
            }
            seq_names.emplace_back(seq_name);
        }
        return seq_names;
    }

    std::vector<hts_pos_t> get_seq_lengths(const faidx_t* faidx)
    {
        std::vector<hts_pos_t> seq_lengths;
        const auto size = faidx_nseq(faidx);
        seq_lengths.reserve(size);

        for (int i = 0; i < size; i++) {
            seq_lengths.emplace_back(faidx_seq_len64(faidx, faidx_iseq(faidx, i)));
        }
        return seq_lengths;
    }

    faidx_t* get_faidx(const std::string& file_name)
    {
        if (!exists(boost::filesystem::path(std::string(fai_path(file_name.c_str()))))) {
            BOOST_LOG_TRIVIAL(fatal) << "FAI not found!";
            abort_mpi();
        }
        BOOST_LOG_TRIVIAL(debug) << "Loading existing FAI...";
        auto* const faidx = fai_load_format(file_name.c_str(), FAI_FASTA);
        if (faidx == nullptr) {
            BOOST_LOG_TRIVIAL(fatal) << "Loading FAI failed!";
            abort_mpi();
        }
        return faidx;
    }

} // namespace

FaidxFetch::~FaidxFetch() { fai_destroy(faidx_); }

FaidxFetch::FaidxFetch(const std::string& file_name)
    : FaidxFetch(get_faidx(file_name))
{
}

FaidxFetch::FaidxFetch(faidx_t* faidx)
    : BaseFastaFetch(get_seq_names(faidx), get_seq_lengths(faidx))
    , faidx_(faidx)
{
}
std::string FaidxFetch::fetch(const size_t seq_id, const hts_pos_t start, const hts_pos_t end)
{
#ifdef CEU_CM_IS_DEBUG
    if (seq_id >= seq_names_.size()) {
        BOOST_LOG_TRIVIAL(fatal) << "InMemoryFastaFetch::fetch: Requested seq_id " << seq_id
                                 << " is out of bounds for total sequences of " << seq_names_.size() << ".";
        abort_mpi();
    }
    const hts_pos_t max_len = faidx_seq_len64(faidx_, seq_names_[seq_id].c_str());

    if (end > max_len || start > max_len || start < 0 || end < 0 || start > end) {
        BOOST_LOG_TRIVIAL(fatal) << "InMemoryFastaFetch::fetch: Requested range [" << start << ", " << end
                                 << ") is out of bounds for sequence of length " << max_len << ".";
        abort_mpi();
    }
#endif
    hts_pos_t len = 0;
    auto* const cfetch_str = faidx_fetch_seq64(faidx_, seq_names_[seq_id].c_str(), start, end - 1, &len);
    if (cfetch_str == nullptr || len != end - start) {
        BOOST_LOG_TRIVIAL(fatal) << "FaidxFetch failed at " << seq_names_[seq_id].c_str() << ":" << start << "-" << end
                                 << "!";
        abort_mpi();
    }
    const auto rets = std::string(cfetch_str, len);
    std::free(cfetch_str);
    return rets;
}
} // namespace labw::art_modern
