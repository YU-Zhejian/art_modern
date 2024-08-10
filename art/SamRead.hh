#pragma once

#include <string>
namespace labw {
namespace art_modern {

    class SamRead {
    public:
        std::string qname;
        int flag = 0;
        std::string rname;
        long pos = 0;
        int mapQ = 99;
        std::string cigar;
        std::string rNext = "*";
        long pNext = 0;
        int tLen = 0;
        std::string seq;
        std::string qual;
        void reverse_comp();
        void printRead(std::ostream& fout) const;
    };

} // namespace art_modern
} // namespace labw