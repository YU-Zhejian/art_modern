#pragma once

#include <htslib/faidx.h>
#include <mutex>

#include "BaseFastaFetch.hh"

namespace labw {
namespace art_modern {

    class FaidxFetch : public BaseFastaFetch {
    public:
        explicit FaidxFetch(faidx_t* faidx);
        explicit FaidxFetch(const std::string& file_name);
        std::string fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end) override;
        ~FaidxFetch() override;

    private:
        faidx_t* faidx_;
        char* cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end);
        std::mutex mutex_;
    };
}
}
