#include "FastaStreamBatcher.hh"
#include <boost/log/trivial.hpp>

namespace labw {
namespace art_modern {
    FastaStreamBatcher::FastaStreamBatcher(std::size_t batch_size, std::istream& stream)
        : batch_size_(batch_size)
        , fasta_iterator_(stream)
    {
    }
    InMemoryFastaFetch FastaStreamBatcher::fetch()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        std::unordered_map<std::string, std::string> fasta_map;
        std::string fetch_s;
        std::string fetch_e;
        while (fasta_map.size() < batch_size_) {
            try {
                auto fasta_record = fasta_iterator_.next();
                if (fetch_s.empty()) {
                    fetch_s = fasta_record.id;
                }
                fetch_e = fasta_record.id;
                fasta_map.emplace(fasta_record.id, fasta_record.sequence);
            } catch (EOFException&) {
                break;
            }
        }
        BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << fetch_s << " to " << fetch_e << " (" << fasta_map.size()
                                << "） created";
        return InMemoryFastaFetch(fasta_map);
    }

    InMemoryFastaFetch InMemoryFastaStreamBatcher::fetch()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        std::unordered_map<std::string, std::string> fasta_map;
        std::string fetch_s;
        std::string fetch_e;
        while (current_index_ < stream_->num_seqs() && fasta_map.size() < batch_size_) {
            auto seq_name = stream_->seq_names_[current_index_];
            fasta_map.emplace(seq_name, stream_->fetch(seq_name, 0, stream_->seq_len(seq_name)));
            if (fetch_s.empty()) {
                fetch_s = seq_name;
            }
            fetch_e = seq_name;
            current_index_++;
        }
        BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << fetch_s << " to " << fetch_e << " (" << fasta_map.size()
                                << "） created";
        return InMemoryFastaFetch(fasta_map);
    }

    InMemoryFastaStreamBatcher::InMemoryFastaStreamBatcher(std::size_t batch_size, BaseFastaFetch* stream)
        : batch_size_(batch_size)
        , current_index_(0)
        , stream_(stream)
    {
    }
}
}