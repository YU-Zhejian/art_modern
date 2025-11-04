#include "libam_support/ref/walker/GenomeWalker.hh"

namespace labw::art_modern {
void GenomeWalker::attach(const std::size_t seq_id, const hts_pos_t start, const bool on_pos_strand) {
    seq_id_ = seq_id;
    position_ = start;
    on_pos_strand_ = on_pos_strand;
}
std::size_t GenomeWalker::get_seq_id() const {
    return seq_id_;
}
hts_pos_t GenomeWalker::get_position() const {
    return position_;
}
bool GenomeWalker::is_on_pos_strand() const {
    return on_pos_strand_;
}
void GenomeWalker::clear() {
    seq_id_ = 0;
    position_ = 0;
    on_pos_strand_ = true;
}
void GenomeWalker::forward(const hts_pos_t offset) {
    if (on_pos_strand_) {
        position_ += offset;
    } else {
        position_ -= offset;
    }
}
bool GenomeWalker::reached_genome_boundary() const {
    if (on_pos_strand_) {
        return position_ >= fasta_fetcher_->seq_len(seq_id_);
    } else {
        return position_ <= 0;
    }
}
GenomeWalker::GenomeWalker(std::shared_ptr<BaseFastaFetch>&& fasta_fetcher)
    : fasta_fetcher_(std::move(fasta_fetcher)) {}
GenomeWalker::GenomeWalker(const std::shared_ptr<BaseFastaFetch>& fasta_fetcher)
    : fasta_fetcher_(fasta_fetcher) {}


} // namespace labw::art_modern
