#include "gffread_bed_utils.hh"

#include "TidNotFound.hh"

#include <htslib/sam.h>

#include <istream>
#include <string>
#include <unordered_map>
namespace labw::art_modern {
std::unordered_map<std::string, Transcript> read_gffutils_bed(sam_hdr_t* thdr, sam_hdr_t* ghdr, std::istream& in)
{
    std::unordered_map<std::string, Transcript> transcript_id_to_transcript_map;
    std::string line;

    while (std::getline(in, line)) {
        try {
            auto transcript = Transcript::from_gffread_bed_line(line, thdr, ghdr);
            transcript_id_to_transcript_map.emplace(transcript.transcript_id, transcript);
        } catch (TidNotFound& /** ignored **/) {
            // ignored
        }
    }
    return transcript_id_to_transcript_map;
}
} // namespace labw::art_modern
