#pragma once

#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstddef>
#include <string>

#include <htslib/faidx.h>
#include <htslib/hts.h>

namespace labw::art_modern {

/**
 * Please note that this method have no lock.
 *
 * That is, all threads MUST own their own FaidxFetch.
 */
class FaidxFetch : public BaseFastaFetch {
public:
    DELETE_COPY(FaidxFetch)
    DELETE_MOVE(FaidxFetch)

    explicit FaidxFetch(faidx_t* faidx);
    explicit FaidxFetch(const std::string& file_name);
    std::string fetch(size_t seq_id, hts_pos_t start, hts_pos_t end) override;
    ~FaidxFetch() override;

private:
    faidx_t* faidx_;
    char* cfetch_(const char* seq_name, hts_pos_t start, hts_pos_t end) const;
};
} // namespace labw::art_modern
