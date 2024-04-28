//
// Created by yuzj on 3/9/24.
//

#include "fasta_parser.hh"
#define BOOST_TEST_MODULE test_fasta_parser
#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace std;
using namespace labw::art_modern;

BOOST_AUTO_TEST_CASE(test_fasta_parser_1)
{
    auto iss = std::istringstream(">chr1\r\nAAAA\nCCC\n>chr2\r\n>chr3\nTTTT\n\n");
    auto fai = FastaIterator(iss);
    vector<string> chrNames = { "chr1", "chr2", "chr3" };
    vector<string> seqs = { "AAAACCC", "", "TTTT" };
    int i = 0;
    while (true) {
        FastaRecord fa_record;
        try {
            fa_record = fai.next();
            BOOST_TEST(fa_record.sequence == seqs[i]);
            BOOST_TEST(fa_record.id == chrNames[i]);
        } catch (EOFException&) {
            break;
        }
        i++;
    }
}
