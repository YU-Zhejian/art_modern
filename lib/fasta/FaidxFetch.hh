//
// Created by yuzj on 24-8-11.
//

#ifndef ART_MODERN_LIB_FASTA_FAIDXFETCH_HH
#define ART_MODERN_LIB_FASTA_FAIDXFETCH_HH

#include <mutex>
#include <htslib/faidx.h>

#include "FastaFetch.hh"

namespace labw{
namespace art_modern {

class FaidxFetch : public FastaFetch {
public:
  explicit FaidxFetch(const std::string& file_name);
  std::string fetch(const std::string &seq_name, hts_pos_t start, hts_pos_t end) override;
  char *cfetch(const char *seq_name, hts_pos_t start, hts_pos_t end) override;
  ~FaidxFetch() override;
private:
  faidx_t* faidx_;
  std::mutex mutex_;
};
}}

#endif //ART_MODERN_LIB_FASTA_FAIDXFETCH_HH
