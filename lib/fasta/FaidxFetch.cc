#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>

#include "FaidxFetch.hh"
#include "utils/mpi_utils.hh"

namespace labw {
namespace art_modern {

    char* FaidxFetch::cfetch_(const char* seq_name, hts_pos_t start, hts_pos_t end)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto reg = boost::format("%s:%d-%d") % seq_name % (start + 1) % end;
        hts_pos_t pos;
        auto rets = fai_fetch64(faidx_, reg.str().c_str(), &pos);
        if (!rets) {
            BOOST_LOG_TRIVIAL(fatal) << "FaidxFetch failed at " << seq_name << ":" << start << "-" << end << "!";
            abort_mpi();
        }
        return rets;
    }
    FaidxFetch::~FaidxFetch() { fai_destroy(faidx_); }

    std::tuple<std::vector<std::string>, std::vector<hts_pos_t>> get_seq_names_lengths(const faidx_t* faidx)
    {
        std::vector<std::string> seq_names;
        std::vector<hts_pos_t> seq_lengths;
        auto size = faidx_nseq(faidx);
        seq_names.reserve(size);
        seq_lengths.reserve(size);

        for (int i = 0; i < size; i++) {
            auto seq_name = faidx_iseq(faidx, i);
            if (!seq_name) {
                BOOST_LOG_TRIVIAL(fatal) << "Sequence name of seq " << i << " is null!";
                abort_mpi();
            }
            seq_names.emplace_back(seq_name);
            seq_lengths.emplace_back(faidx_seq_len(faidx, seq_name));
        }
        return std::tie(seq_names, seq_lengths);
    }

    faidx_t* get_faidx(const std::string& file_name)
    {
        auto seq_file_fai_path = std::string(fai_path(file_name.c_str()));
        if (!boost::filesystem::exists(boost::filesystem::path(seq_file_fai_path))) {
            BOOST_LOG_TRIVIAL(fatal) << "FAI not found!";
            abort_mpi();
        } else {
            BOOST_LOG_TRIVIAL(info) << "Loading existing FAI...";
            auto faidx = fai_load_format(file_name.c_str(), FAI_FASTA);
            if (!faidx) {
                BOOST_LOG_TRIVIAL(fatal) << "Loading FAI failed!";
                abort_mpi();
            }
            return faidx;
        }
    }

    FaidxFetch::FaidxFetch(const std::string& file_name)
        : FaidxFetch(get_faidx(file_name))
    {
    }

    FaidxFetch::FaidxFetch(faidx_t* faidx)
        : BaseFastaFetch(get_seq_names_lengths(faidx))
        , faidx_(faidx)
    {
    }
    std::string FaidxFetch::fetch(size_t seq_id, hts_pos_t start, hts_pos_t end)
    {
        auto cfetch_str = cfetch_(seq_names_[seq_id].c_str(), start, end);
        auto rets = std::string(cfetch_str);
        std::free(cfetch_str);
        return rets;
    }
}
}