/**
 * Copyright 2026 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include "libam_support/bam/BamOptions.hh"

#include "libam_support/utils/seq_utils.hh"

#include <boost/log/trivial.hpp>

#include <string>
#include <vector>

namespace labw::art_modern {

void BamOptions::log_(const std::string& name) const
{
    std::vector<std::string> tags = {};
    if (with_tag_MD) {
        tags.emplace_back("MD");
    }
    if (with_tag_NM) {
        tags.emplace_back("NM");
    }
    if (with_tag_OA) {
        tags.emplace_back("OA");
    }
    BOOST_LOG_TRIVIAL(info) << name << ": SAM/BAM Output Options: use_m=" << use_m << ", write_bam=" << write_bam
                            << ", hts_io_threads=" << hts_io_threads << ", compress_level=" << compress_level
                            << ", tags=[" << join(tags, ",") << "], no_qual=" << no_qual << ".";
}
} // namespace labw::art_modern
