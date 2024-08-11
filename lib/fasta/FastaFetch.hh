#ifndef ART_MODERN_LIB_FASTAFETCH_HH
#define ART_MODERN_LIB_FASTAFETCH_HH
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <map>
#include <string>

namespace labw {
namespace art_modern {
    class FastaFetch {
    public:
        virtual std::string fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end);
        /**
         *
         * @param seq_name
         * @param start
         * @param end
         * @return The returned string should be manually freed.
         */
        virtual char* cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end);
        virtual ~FastaFetch();
        void update_sam_header(sam_hdr_t* header);
        hts_pos_t seq_len(const std::string& seq_name);
        size_t num_seqs();

        std::map<std::string, hts_pos_t, std::less<>> seq_lengths_;
    };
}
}

#endif // ART_MODERN_LIB_FASTAFETCH_HH
