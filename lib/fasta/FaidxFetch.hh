#pragma once

#include <htslib/faidx.h>
#include <mutex>

#include "BaseFastaFetch.hh"

namespace labw::art_modern {

/**
 * Please note that this method have no lock.
 *
 * That is, all threads MUST own their own FaidxFetch.
 */
    class FaidxFetch : public BaseFastaFetch {
    public:
        explicit FaidxFetch(faidx_t* faidx);
        explicit FaidxFetch(const std::string& file_name);
        std::string fetch(size_t seq_id, hts_pos_t start, hts_pos_t end) override;
        ~FaidxFetch() override;

        FaidxFetch(const FaidxFetch&) = delete;
        FaidxFetch(FaidxFetch&&) = delete;
        FaidxFetch& operator=(const FaidxFetch&) = delete;
        FaidxFetch& operator=(FaidxFetch&&) = delete;

    private:
        faidx_t* faidx_;
        char* cfetch_(const char* seq_name, hts_pos_t start, hts_pos_t end);
    };
}

