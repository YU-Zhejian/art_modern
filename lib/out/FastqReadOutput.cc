#include "FastqReadOutput.hh"
namespace labw {
namespace art_modern {

    void FastqReadOutput::writeSE(const PairwiseAlignment& pwa)
    {
        stream_.write("@");
        stream_.write(pwa.read_name);
        stream_.write("\n");
        stream_.write(pwa.query);
        stream_.write("+\n");
        stream_.write(pwa.qual);
        stream_.write("\n");
    }

    void FastqReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        stream_.write("@");
        stream_.write(pwa1.read_name);
        stream_.write("/1\n");
        stream_.write(pwa1.query);
        stream_.write("+\n");
        stream_.write(pwa1.qual);
        stream_.write("\n@");
        stream_.write(pwa2.read_name);
        stream_.write("/2\n");
        stream_.write(pwa2.query);
        stream_.write("+\n");
        stream_.write(pwa2.qual);
        stream_.write("\n");
    }

    FastqReadOutput::~FastqReadOutput() = default;
    FastqReadOutput::FastqReadOutput(ThreadSafeFileStream& stream)
        : stream_(stream)
    {
    }
}
}