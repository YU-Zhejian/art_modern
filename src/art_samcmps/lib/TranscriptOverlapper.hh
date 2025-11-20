#pragma once

#include "art_samcmps/lib/AlignmentMiscInfo.hh"

#include "libam_support/utils/class_macros_utils.hh"

#include <IITree.h>

#include <htslib/hts.h>

#include <cstdlib>
#include <istream>
#include <string>
#include <vector>

namespace labw::art_modern {
class TranscriptOverlapper {

    std::vector<std::string> transcript_ids;
    std::vector<std::string> transcript_contig_names;
    std::vector<AlignmentMiscInfo> amis;
    IITree<hts_pos_t, std::size_t> transcript_ranges;

public:
    DELETE_COPY(TranscriptOverlapper)

    DELETE_MOVE(TranscriptOverlapper)

    ~TranscriptOverlapper() = default;

    explicit TranscriptOverlapper(std::istream& bed_stream);

    std::vector<AlignmentMiscInfo*> get_overlapping_transcripts(const std::string& contig, int start, int end);
};
} // namespace labw::art_modern
