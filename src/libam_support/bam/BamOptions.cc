#include "libam_support/bam/BamOptions.hh"

#include "libam_support/utils/seq_utils.hh"

#include <boost/log/trivial.hpp>

#include <string>
#include <vector>

namespace labw::art_modern {

void BamOptions::log_(const std::string& name) const
{
    std::vector<std::string> tags = {};
    if( with_tag_MD) {
        tags.emplace_back( "MD" );
    }
    if( with_tag_NM) {
        tags.emplace_back( "NM" );
    }
        if( with_tag_OA) {
                tags.emplace_back( "OA" );
        }
    BOOST_LOG_TRIVIAL( info ) << name << ": SAM/BAM Output Options: use_m=" << use_m << ", write_bam=" << write_bam
                      << ", hts_io_threads=" << hts_io_threads << ", compress_level=" << compress_level
                      << ", tags=[" << join( tags, "," ) << "], no_qual=" << no_qual << ".";
}
} // namespace labw::art_modern
