//
// Created by yuzj on 24-8-11.
//

#include <utility>

#include "InMemoryFastaFetch.hh"


namespace labw {
namespace art_modern {
InMemoryFastaFetch::InMemoryFastaFetch(std::map<std::string, std::string> seq_map_): seq_map_(std::move(seq_map_)) {
  for(auto const & pair: seq_map_ ){
    seq_lengths_.emplace(pair.first, pair.second.size());
  }
}
std::string InMemoryFastaFetch::fetch(const std::string &seq_name, hts_pos_t start, hts_pos_t end) {
  return seq_map_.at(seq_name).substr(start, end);
}
char *InMemoryFastaFetch::cfetch(const char *seq_name, hts_pos_t start, hts_pos_t end) {
  auto fetch_str = fetch(seq_name, start, end);
  auto rets = (char*) calloc(fetch_str.size() + 1, sizeof (char));
  strncpy(rets, fetch_str.c_str(), fetch_str.size());
  return rets;
}

InMemoryFastaFetch::~InMemoryFastaFetch() =default;
}}