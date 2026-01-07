// Using NCBI NGS C++ API to read a local PE SRA file

#include <ngs/Read.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>

#include <ncbi-vdb/NGS.hpp>

#include <cstdlib>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_sra_file>" << std::endl;
        return EXIT_FAILURE;
    }
    const char* sra_path = argv[1];
    try {
        ngs::ReadCollection run = ncbi::NGS::openReadCollection(sra_path);
        auto const nreads = run.getReadCount();
        std::cout << "Total Reads: " << nreads << std::endl;

        auto reads = run.getReadRange(1, nreads);

        while (reads.nextRead()) {
            std::cout << "Read ID: " << reads.getReadName() << std::endl;
            reads.nextFragment();
            std::cout << "Sequence 1: " << reads.getFragmentBases().toString() << std::endl;
            std::cout << "Quality 1: " << reads.getFragmentQualities().toString() << std::endl;
            reads.nextFragment();
            std::cout << "Sequence 2: " << reads.getFragmentBases().toString() << std::endl;
            std::cout << "Quality 2: " << reads.getFragmentQualities().toString() << std::endl;
        }
    } catch (const ngs::ErrorMsg& e) {
        std::cerr << "NCBI NGS Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (const std::exception& e) {
        std::cerr << "Standard Exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
