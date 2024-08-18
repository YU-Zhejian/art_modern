#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>

#include "FaidxFetch.hh"
namespace labw {
namespace art_modern {

    std::string FaidxFetch::fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end)
    {
        auto cfetch_str = cfetch(seq_name.c_str(), start, end);
        auto rets = std::string(cfetch_str);
        free(cfetch_str);
        return rets;
    }

    char* FaidxFetch::cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto reg = boost::format("%s:%d-%d") % seq_name % (start + 1) % end;
        auto pos = (hts_pos_t*)malloc(sizeof(hts_pos_t));
        auto rets = fai_fetch64(faidx_, reg.str().c_str(), pos);
        if (!rets) {
            BOOST_LOG_TRIVIAL(fatal) << "FaidxFetch failed at " << seq_name << ":" << start << "-" << end << "!";
            exit(EXIT_FAILURE);
        }
        free(pos);
        return rets;
    }
    FaidxFetch::~FaidxFetch() { fai_destroy(faidx_); }
    FaidxFetch::FaidxFetch(const std::string& file_name)
    {
        auto seq_file_fai_path = std::string(fai_path(file_name.c_str()));
        if (!boost::filesystem::exists(boost::filesystem::path(seq_file_fai_path))) {
            BOOST_LOG_TRIVIAL(fatal) << "FAI not found!";
            exit(EXIT_FAILURE);
        } else {
            BOOST_LOG_TRIVIAL(info) << "Loading existing FAI...";
            faidx_ = fai_load_format(file_name.c_str(), FAI_FASTA);
            if (!faidx_) {
                BOOST_LOG_TRIVIAL(fatal) << "Loading FAI failed!";
                exit(EXIT_FAILURE);
            }
        }
        for (int i = 0; i < faidx_nseq(faidx_); i++) {
            auto seq_name = faidx_iseq(faidx_, i);
            seq_lengths_.emplace(std::string(seq_name), faidx_seq_len(faidx_, seq_name));
        }
    }
}
}