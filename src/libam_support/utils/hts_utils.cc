#include "libam_support/utils/hts_utils.hh"

#include <htslib/bgzf.h> // NOLINT
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <cstddef>

namespace labw::art_modern {

std::size_t hts_tell(htsFile* fp)
{
    if (fp == nullptr) {
        return 0;
    }
    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        return htell(fp->fp.bgzf->fp);
    case cram:
        return -1; // TODO: no support for cram
    case empty_format:
    case text_format:
    case bed:
    case fasta_format:
    case fastq_format:
    case sam:
    case vcf:
        if (fp->format.compression != no_compression) {
            return htell(fp->fp.bgzf->fp);
        } else {
            return htell(fp->fp.hfile);
        }
    default:
        break;
    }
    return 0;
}

} // namespace labw::art_modern