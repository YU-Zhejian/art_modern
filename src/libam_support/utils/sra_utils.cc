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
#include "art_modern_config.h" // NOLINT: for various flags

#include "libam_support/utils/sra_utils.hh"

#include "libam_support/utils/mpi_utils.hh"

#ifdef WITH_NCBI_NGS
#include <ncbi-vdb/NGS.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/itf/ErrorMsg.hpp>
#endif

#include <boost/log/trivial.hpp>

#include <cstring>
#include <fstream>
#include <ios>
#include <string>

namespace labw::art_modern {
bool detect_sra([[maybe_unused]] const std::string& file_path)
{
#ifdef WITH_NCBI_NGS
    auto ifs = std::ifstream(file_path, std::ios::in | std::ios::binary);
    char signature[8] = {};
    ifs.read(signature, 8);
    ifs.seekg(0);
    ifs.close();
    return std::memcmp(signature, SRA_SIGNATURE, 8) == 0;
#else
    BOOST_LOG_TRIVIAL(error) << "NCBI NGS support is not enabled.";
    abort_mpi();
#endif
}

void assert_is_sra([[maybe_unused]] const std::string& file_path)
{
#ifdef WITH_NCBI_NGS
    // Assert that NGS can open the file
    try {
        ngs::ReadCollection const rcngs = ncbi::NGS::openReadCollection(file_path);
    } catch (const ngs::ErrorMsg& e) {
        BOOST_LOG_TRIVIAL(error) << "NGS error when opening SRA file: " << e.what();
        abort_mpi();
    }
#else
    BOOST_LOG_TRIVIAL(error) << "NCBI NGS support is not enabled.";
    abort_mpi();
#endif
}
} // namespace  labw::art_modern
