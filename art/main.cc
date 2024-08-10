#include <boost/asio.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>

#include <fstream>
#include <string>

#include <htslib/faidx.h>

#include "ArtContig.hh"
#include "Empdist.hh"
#include "art_modern_constants.hh"
#include "fasta_parser.hh"
#include "main_fn.hh"

using namespace std;
using namespace labw::art_modern;

int main(int argc, char* argv[])
{
    vector<string> args;
    for (auto i = 0; i < argc; i++) {
        args.emplace_back(argv[i]);
    }

    print_banner();
    ArtParams art_params;
    art_params.parse_args(args);
    art_params.validate_args();
    art_params.print_params();

    auto qdist = art_params.read_emp();
    BOOST_LOG_TRIVIAL(info) << "Checking FAI...";
    auto seq_file_fai_path = fai_path(art_params.seq_file.c_str());
    if (!boost::filesystem::exists(boost::filesystem::path(seq_file_fai_path))) {
        if (!fai_build(art_params.seq_file.c_str())) {
            return 1;
        }
    }
    free(seq_file_fai_path);

    auto fa_reader = std::ifstream(art_params.seq_file);
    FastaIterator fai(fa_reader);
    int num_cores = 1;
    if (art_params.parallel == PARALLEL_ALL) {
        num_cores = static_cast<int>(boost::thread::hardware_concurrency());
    } else if (art_params.parallel != PARALLEL_DISABLE) {
        num_cores = art_params.parallel;
    }
    BOOST_LOG_TRIVIAL(info) << "Using " << num_cores << " cores.";
    boost::asio::thread_pool pool(num_cores);
    while (true) {
        string id;
        string ref_seq;
        try {
            auto fasta_record = fai.next();
            id = fasta_record.id;
            ref_seq = fasta_record.sequence;
        } catch (EOFException&) {
            break;
        }
        auto func = [art_params, id, ref_seq, qdist] {
            double sequencing_depth;
            if (art_params.uniform_sequencing_depth != 0.0) {
                sequencing_depth = art_params.uniform_sequencing_depth;
            } else {
                auto find_result = art_params.sequencing_depth.find(id);
                if (find_result == art_params.sequencing_depth.end()) {
                    // BOOST_LOG_TRIVIAL(warning) << "Warning: depth info for contig '" << id << "' not found!";
                    sequencing_depth = 0;
                } else {
                    sequencing_depth = find_result->second;
                }
            }
            generate_all(id, ref_seq, art_params, qdist, sequencing_depth);
        };
        if (art_params.parallel_on_read || art_params.parallel == PARALLEL_DISABLE) {
            func();
        } else {
            post(pool, func);
        }
    }
    pool.join();
    return EXIT_SUCCESS;
}
