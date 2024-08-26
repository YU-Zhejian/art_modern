#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>

#include "FaidxFetch.hh"
namespace labw {
namespace art_modern {

    std::string FaidxFetch::fetch(const std::string& seq_name, const hts_pos_t start, const hts_pos_t end)
    {
        auto cfetch_str = cfetch_(seq_name.c_str(), start, end);
        auto rets = std::string(cfetch_str);
        free(cfetch_str);
        return rets;
    }

    char* FaidxFetch::cfetch_(const char* seq_name, hts_pos_t start, hts_pos_t end)
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

    std::unordered_map<std::string, hts_pos_t> get_seq_lengths(const faidx_t* faidx)
    {
        std::unordered_map<std::string, hts_pos_t> seq_lengths;
        for (int i = 0; i < faidx_nseq(faidx); i++) {
            auto seq_name = faidx_iseq(faidx, i);
            seq_lengths.emplace(std::string(seq_name), faidx_seq_len(faidx, seq_name));
        }
        return seq_lengths;
    }

    faidx_t* get_faidx(const std::string& file_name)
    {
        auto seq_file_fai_path = std::string(fai_path(file_name.c_str()));
        if (!boost::filesystem::exists(boost::filesystem::path(seq_file_fai_path))) {
            BOOST_LOG_TRIVIAL(fatal) << "FAI not found!";
            exit(EXIT_FAILURE);
        } else {
            BOOST_LOG_TRIVIAL(info) << "Loading existing FAI...";
            auto faidx = fai_load_format(file_name.c_str(), FAI_FASTA);
            if (!faidx) {
                BOOST_LOG_TRIVIAL(fatal) << "Loading FAI failed!";
                exit(EXIT_FAILURE);
            }
            return faidx;
        }
    }

    FaidxFetch::FaidxFetch(const std::string& file_name)
        : FaidxFetch(get_faidx(file_name))
    {
    }

    FaidxFetch::FaidxFetch(faidx_t* faidx)
        : BaseFastaFetch(get_seq_lengths(faidx))
        , faidx_(faidx)
    {
    }
}
}