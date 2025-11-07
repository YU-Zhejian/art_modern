#pragma once

#include "Transcript.hh"

#include <htslib/hts.h>

#include <istream>
#include <string>
#include <unordered_map>
namespace labw::art_modern {
std::unordered_map<std::string, Transcript> read_gffutils_bed(sam_hdr_t* thdr, sam_hdr_t* ghdr, std::istream& in);
} // namespace labw::art_modern
