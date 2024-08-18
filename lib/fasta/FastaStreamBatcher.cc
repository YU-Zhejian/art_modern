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
}
}