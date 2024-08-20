#include <boost/algorithm/string/join.hpp>
#include <boost/asio.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>
#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <fstream>
#include <string>

#include <htslib/faidx.h>

#include "ArtContig.hh"
#include "CExceptionsProxy.hh"
#include "art_modern_constants.hh"
#include "fasta/InMemoryFastaFetch.hh"
#include "fasta/fasta_parser.hh"
#include "global_variables.hh"
#include "main_fn.hh"
#include "out/BamReadOutput.hh"
#include "out/FastqReadOutput.hh"

using namespace std;
using namespace labw::art_modern;

int main(int argc, char* argv[])
{
#ifdef WITH_BOOST_TIMER
    boost::timer::auto_cpu_timer t;
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif
    for (auto i = 0; i < argc; i++) {
        args.emplace_back(argv[i]);
    }

    print_banner();
    ArtParams art_params;
    art_params.parse_args(args);
    art_params.validate_args();
    art_params.print_params();

    art_params.read_emp();
    if (!art_params.stream) {
        BOOST_LOG_TRIVIAL(info) << "HTSLib parser requested. Checking FAI...";
        auto seq_file_fai_path = string(fai_path(art_params.seq_file.c_str()));
        if (!boost::filesystem::exists(boost::filesystem::path(seq_file_fai_path))) {
            BOOST_LOG_TRIVIAL(info) << "Building missing FAI...";
            CExceptionsProxy::requires_numeric(
                fai_build(art_params.seq_file.c_str()), USED_HTSLIB_NAME, "Failed to build FAI");
        } else {
            BOOST_LOG_TRIVIAL(info) << "Loading existing FAI...";
            CExceptionsProxy::requires_not_null(
                fai_load_format(art_params.seq_file.c_str(), FAI_FASTA), USED_HTSLIB_NAME, "Failed to load FAI");
        }
    }

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
        string contig_name;
        string ref_seq;
        try {
            auto fasta_record = fai.next();
            contig_name = fasta_record.id;
            ref_seq = fasta_record.sequence;
        } catch (EOFException&) {
            break;
        }

        auto output_dispatcher = art_params.get_output_dispatcher();
        auto fasta_fetch = std::make_shared<InMemoryFastaFetch>(contig_name, ref_seq);
        SamReadOutputOptions sam_opts;
        sam_opts.write_bam = false;
        sam_opts.PG_CL = boost::algorithm::join(args, " ");

        auto func = [art_params, contig_name, ref_seq, qdist, &output_dispatcher] {
            double sequencing_depth;
            if (art_params.uniform_sequencing_depth != 0.0) {
                sequencing_depth = art_params.uniform_sequencing_depth;
            } else {
                auto find_result = art_params.sequencing_depth.find(contig_name);
                if (find_result == art_params.sequencing_depth.end()) {
                    // BOOST_LOG_TRIVIAL(warning) << "Warning: depth info for contig '" <<
                    // id << "' not found!";
                    sequencing_depth = 0;
                } else {
                    sequencing_depth = find_result->second;
                }
            }
            generate_all(contig_name, ref_seq, art_params, qdist, sequencing_depth, output_dispatcher);
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
