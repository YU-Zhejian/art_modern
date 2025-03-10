#include "test_adaptor.h"

#define BOOST_TEST_MODULE test_fasta_parser // NOLINT

#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/ref/fetch/FaidxFetch.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/ref/parser/fasta_parser.hh"

#include <boost/test/unit_test.hpp>

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace labw::art_modern;

namespace {

void test_fasta(std::unique_ptr<BaseFastaFetch> fastaFetch)
{
    BOOST_TEST(fastaFetch->num_seqs() == 5);
    BOOST_TEST(fastaFetch->fetch(2, 2, 15) == "TANNTGNATNATG");
    BOOST_TEST(fastaFetch->fetch(2, 2, 16) == "TANNTGNATNATGN");
    BOOST_TEST(fastaFetch->fetch(1, 0, fastaFetch->seq_len(1))
        == "NNNNNNNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTTCTNNNNNNNNN");
    BOOST_TEST(fastaFetch->fetch(4, 0, fastaFetch->seq_len(4)) == "CTA");
    BOOST_TEST(fastaFetch->fetch(3, 0, fastaFetch->seq_len(3)) == "AAAAAAAAAACCCCCC");
    BOOST_TEST(fastaFetch->fetch(0, 0, 1) == "N");
    BOOST_TEST(fastaFetch->fetch(0, 26, 29) == "CCA");
    BOOST_TEST(fastaFetch->fetch(0, 28, 29) == "A");
    BOOST_TEST(fastaFetch->fetch(0, 5, 29) == "NNNNNNNNNNATCGTTACGTACCA");
    BOOST_TEST(fastaFetch->fetch(0, 5, 63) == "NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT");
    BOOST_TEST(fastaFetch->seq_len(0) == 154);
    BOOST_TEST(fastaFetch->seq_len(1) == 74);
}

} // namespace

BOOST_AUTO_TEST_CASE(test_fasta_parser_1)
{
    std::istringstream iss(">chr1\r\nAAAA\nCCC\n>chr2\r\n>chr3\nTTTT\n\n");
    FastaIterator fai(iss);
    std::vector<std::string> chrNames = { "chr1", "chr2", "chr3" };
    std::vector<std::string> seqs = { "AAAACCC", "", "TTTT" };
    int i = 0;
    while (true) {
        try {
            const auto& [id, sequence] = fai.next();
            BOOST_TEST(id == chrNames[i]);
            BOOST_TEST(sequence == seqs[i]);
        } catch (EOFException&) {
            break;
        }
        i++;
    }
    BOOST_TEST(i == chrNames.size());
}

BOOST_AUTO_TEST_CASE(test_faidx_fetch)
{
    auto faidx_fetch = std::make_unique<FaidxFetch>(TEST_RESOURCES_PATH "test.fasta");
    test_fasta(std::move(faidx_fetch));
}

BOOST_AUTO_TEST_CASE(test_in_memory_fetch)
{
    auto in_memory_fasta_fetch = std::make_unique<InMemoryFastaFetch>(TEST_RESOURCES_PATH "test.fasta");
    test_fasta(std::move(in_memory_fasta_fetch));
}
