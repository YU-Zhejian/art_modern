#include <utility>

#include "FastaStreamBatcher.hh"

namespace labw {
namespace art_modern {
    FastaStreamBatcher::FastaStreamBatcher(int batch_size, std::istream& stream)
        : batch_size_(batch_size)
        , fasta_iterator_(stream)
    {
    }
    InMemoryFastaFetch FastaStreamBatcher::fetch()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        std::map<std::string, std::string, std::less<>> fasta_map;
        while (fasta_map.size() < batch_size_) {
            FastaRecord fasta_record;
            try {
                fasta_record = fasta_iterator_.next();
            } catch (EOFException&) {
                break;
            }
            fasta_map.emplace(fasta_record.id, fasta_record.sequence);
        }
        return InMemoryFastaFetch(fasta_map);
    }

    InMemoryFastaFetch InMemoryFastaStreamBatcher::fetch()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        std::map<std::string, std::string, std::less<>> fasta_map;
        while (current_index_ < stream_->num_seqs() && fasta_map.size() < batch_size_) {
            auto seq_name = seq_names_[current_index_];
            fasta_map.emplace(seq_name, stream_->fetch(seq_name, 0, stream_->seq_len(seq_name)));
            current_index_++;
        }
        return InMemoryFastaFetch(fasta_map);
    }

    InMemoryFastaStreamBatcher::InMemoryFastaStreamBatcher(int batch_size, std::shared_ptr<BaseFastaFetch> stream)
        : batch_size_(batch_size)
        , stream_(std::move(stream))
        , current_index_(0)
        , seq_names_(stream_->seq_names())
    {
    }
}
}