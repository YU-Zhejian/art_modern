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
        std::vector<std::string> seq_names;
        std::vector<std::string> seqs;
        seq_names.reserve(batch_size_);
        seqs.reserve(batch_size_);

        std::string fetch_s;
        std::string fetch_e;
        while (seq_names.size() < batch_size_) {
            try {
                auto fasta_record = fasta_iterator_.next();
                if (fetch_s.empty()) {
                    fetch_s = fasta_record.id;
                }
                fetch_e = fasta_record.id;
                seq_names.emplace_back(fasta_record.id);
                seqs.emplace_back(fasta_record.sequence);
            } catch (EOFException&) {
                break;
            }
        }
        BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << fetch_s << " to " << fetch_e << " (" << seq_names.size()
                                << "） created";
        return { seq_names, seqs };
    }

    InMemoryFastaFetch InMemoryFastaStreamBatcher::fetch()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        std::vector<std::string> seq_names;
        std::vector<std::string> seqs;
        seq_names.reserve(batch_size_);
        seqs.reserve(batch_size_);
        while (current_index_ < stream_->num_seqs() && seq_names.size() < batch_size_) {
            auto seq_name = stream_->seq_name(current_index_);
            seq_names.emplace_back(seq_name);
            seqs.emplace_back(stream_->fetch(current_index_, 0, stream_->seq_len(current_index_)));
            current_index_++;
        }
        if (!seq_names.empty()) {
            BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << seq_names[0] << " to " << seq_names[seq_names.size() - 1]
                                    << " (" << seq_names.size() << "） created";
        }

        return { seq_names, seqs };
    }

    InMemoryFastaStreamBatcher::InMemoryFastaStreamBatcher(std::size_t batch_size, BaseFastaFetch* stream)
        : batch_size_(batch_size)
        , current_index_(0)
        , stream_(stream)
    {
    }
}
}