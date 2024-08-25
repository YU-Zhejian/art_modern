//
// Created by yuzj on 3/9/24.
//
#define BOOST_TEST_MODULE test_fasta_parser
#include <boost/test/unit_test.hpp>

#include "fasta/BaseFastaFetch.hh"
#include "fasta/FaidxFetch.hh"
#include "fasta/InMemoryFastaFetch.hh"
#include "fasta/fasta_parser.hh"
#include "test_adaptor.h"

BOOST_AUTO_TEST_CASE(test_fasta_parser_1)
{
    std::istringstream iss(">chr1\r\nAAAA\nCCC\n>chr2\r\n>chr3\nTTTT\n\n");
    labw::art_modern::FastaIterator fai(iss);
    std::vector<std::string> chrNames = { "chr1", "chr2", "chr3" };
    std::vector<std::string> seqs = { "AAAACCC", "", "TTTT" };
    int i = 0;
    while (true) {
        try {
            auto fa_record = fai.next();
            BOOST_TEST(fa_record.id == chrNames[i]);
            BOOST_TEST(fa_record.sequence == seqs[i]);
        } catch (labw::art_modern::EOFException&) {
            break;
        }
        i++;
    }
    BOOST_TEST(i == chrNames.size());
}

void test_fasta(labw::art_modern::BaseFastaFetch* fastaFetch)
{
    BOOST_TEST(fastaFetch->num_seqs() == 5);
    BOOST_TEST(fastaFetch->fetch("chr3", 2, 15) == "TANNTGNATNATG");
    BOOST_TEST(fastaFetch->fetch("chr3", 2, 16) == "TANNTGNATNATGN");
    BOOST_TEST(fastaFetch->fetch("chr2", 0, fastaFetch->seq_len("chr2"))
        == "NNNNNNNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTTCTNNNNNNNNN");
    BOOST_TEST(fastaFetch->fetch("chr6", 0, fastaFetch->seq_len("chr6")) == "CTA");
    BOOST_TEST(fastaFetch->fetch("chr4", 0, fastaFetch->seq_len("chr4")) == "AAAAAAAAAACCCCCC");
    BOOST_TEST(fastaFetch->fetch("chr1", 0, 1) == "N");
    BOOST_TEST(fastaFetch->fetch("chr1", 26, 29) == "CCA");
    BOOST_TEST(fastaFetch->fetch("chr1", 28, 29) == "A");
    BOOST_TEST(fastaFetch->fetch("chr1", 5, 29) == "NNNNNNNNNNATCGTTACGTACCA");
    BOOST_TEST(fastaFetch->fetch("chr1", 5, 63) == "NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT");
}

BOOST_AUTO_TEST_CASE(test_faidx_fetch)
{
    auto faidx_fetch = new labw::art_modern::FaidxFetch(TEST_RESOURCES_PATH "test.fasta");
    test_fasta(faidx_fetch);
    delete faidx_fetch;
}

BOOST_AUTO_TEST_CASE(test_in_memory_fetch)
{
    auto in_memory_fasta_fetch = new labw::art_modern::InMemoryFastaFetch(TEST_RESOURCES_PATH "test.fasta");
    test_fasta(in_memory_fasta_fetch);
    delete in_memory_fasta_fetch;
}
