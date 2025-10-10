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

#include "libam_support/utils/hts_utils.hh"

#include <htslib/bgzf.h> // NOLINT
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h> // NOLINT

#include <cstddef>

namespace labw::art_modern {

std::size_t hts_tell(htsFile* fp)
{
    if (fp == nullptr) {
        return 0;
    }
    switch (fp->format.format) {
    case binary_format: /** fall through */
    case bam: /** fall through */
    case bcf: /** fall through */
        return htell(fp->fp.bgzf->fp);
    case cram: /** fall through */
        return -1;
    case empty_format: /** fall through */
    case text_format: /** fall through */
    case bed: /** fall through */
    case fasta_format: /** fall through */
    case fastq_format: /** fall through */
    case sam: /** fall through */
    case vcf: /** fall through */
    {
        if (fp->format.compression != no_compression) {
            return htell(fp->fp.bgzf->fp);
        }
        return htell(fp->fp.hfile);
    }
    default:
        break;
    }
    return 0;
}
} // namespace labw::art_modern
