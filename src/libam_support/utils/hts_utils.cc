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
    case binary_format: /** fall through */
    case bam: /** fall through */
    case bcf: /** fall through */
        return htell(fp->fp.bgzf->fp);
    case cram: /** fall through */
        return -1; // TODO: no support for cram
    case empty_format: /** fall through */
    case text_format: /** fall through */
    case bed: /** fall through */
    case fasta_format: /** fall through */
    case fastq_format: /** fall through */
    case sam: /** fall through */
    case vcf: /** fall through */
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