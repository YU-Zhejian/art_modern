#include "art_samcmps/lib/TranscriptOverlapper.hh"

#include "art_samcmps/lib/AlignmentMiscInfo.hh"

#include "art_tsam2gsam/lib/TidNotFound.hh"

#include <IITree.h>

#include <cstdlib>
#include <iostream>
#include <istream>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
TranscriptOverlapper::TranscriptOverlapper(std::istream& bed_stream)
{
    std::string line;
    std::size_t line_id = 0;
    while (std::getline(bed_stream, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        try {
            AlignmentMiscInfo ami {};
            ami.from_gffutils_bed_line(line);
            transcript_ids.emplace_back(ami.name);
            transcript_ranges.add(ami.start, ami.end, line_id);
            transcript_contig_names.emplace_back(ami.contig_name);
            line_id++;
            amis.emplace_back(std::move(ami));
        } catch (TidNotFound& e) {
            continue;
        }
    }
    std::cerr << "Loaded " << transcript_ids.size() << " transcripts\n";
    transcript_ranges.index();
}

std::vector<AlignmentMiscInfo*> TranscriptOverlapper::get_overlapping_transcripts(
    const std::string& contig, const int start, const int end)
{
    std::vector<AlignmentMiscInfo*> retl;
    std::vector<std::size_t> id {};

    transcript_ranges.overlap(start, end, id);
    for (const auto& i : id) {
        if (transcript_contig_names[i] != contig) {
            continue;
        }
        retl.emplace_back(&amis[i]);
    }
    return retl;
}
} // namespace labw::art_modern
