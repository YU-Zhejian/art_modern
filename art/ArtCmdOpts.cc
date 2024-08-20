#include "ArtCmdOpts.hh"

#include <boost/log/trivial.hpp>
#include <boost/algorithm/string/join.hpp>
#include "art_modern_constants.hh"
#include "out/BamReadOutput.hh"
#include "out/HeadlessBamReadOutput.hh"
#include "out/FastqReadOutput.hh"
#include "out/OutputDispatcher.hh"
#include "ArtConstants.hh"
#include "ArtVersion.hh"
#include "global_variables.hh"

using namespace std;
namespace po = boost::program_options;

namespace labw
{
namespace art_modern
{

void print_help(const po::options_description& po_desc){
    std::cout << po_desc << endl;
}

po::variables_map generate_vm_while_handling_help_version(const po::options_description& po_desc, int argc, char** argv) noexcept {
    BOOST_LOG_TRIVIAL(info) << "ARGS: " << boost::algorithm::join(args, " ");
    if (argc <= 1) {
        print_help(po_desc);
        exit(EXIT_FAILURE);
    }
    po::variables_map vm_;

    try {
        po::store(po::parse_command_line(static_cast<int>(args.size()), argv, po_desc), vm_);
    } catch (const exception& exp) {
        BOOST_LOG_TRIVIAL(fatal) << exp.what();
        print_help(po_desc);
        exit(EXIT_FAILURE);
    }
    po::notify(vm_);

    if (vm_.count("version")) {
        print_version();
        exit(EXIT_SUCCESS);
    }
    if (vm_.count("help")) {
        print_help(po_desc);
        exit(EXIT_SUCCESS);
    }
    return vm_;
}

ArtParams ArtCmdOpts::parse_args(int argc, char** argv){
        for (int i = 0; i < argc; i++) {
            args.emplace_back(argv[i]);
        }
        auto vm_ = generate_vm_while_handling_help_version(po_desc, argc, argv);

        auto simulation_mode_str = vm_["mode"].as<string>();
    SIMULATION_MODE art_simulation_mode;
        if (simulation_mode_str == "wgs") {
            art_simulation_mode = SIMULATION_MODE::WGS;
        } else if (simulation_mode_str == "trans") {
            art_simulation_mode = SIMULATION_MODE::TRANS;
        } else if (simulation_mode_str == "template") {
            art_simulation_mode = SIMULATION_MODE::TEMPLATE;
        } else {
            BOOST_LOG_TRIVIAL(fatal) << R"(Simulation mode (--mode) should be one of "wgs", "trans" and "template")";
            exit(EXIT_FAILURE);
        }
        auto lib_const_mode_str = vm_["lc"].as<string>();
        ART_LIB_CONST_MODE art_lib_const_mode;
        if (lib_const_mode_str == ART_LIB_CONST_MODE_SE) {
            art_lib_const_mode = ART_LIB_CONST_MODE::SE;
        } else if (lib_const_mode_str == ART_LIB_CONST_MODE_PE) {
            art_lib_const_mode = ART_LIB_CONST_MODE::PE;
        } else if (lib_const_mode_str == ART_LIB_CONST_MODE_MP) {
            art_lib_const_mode = ART_LIB_CONST_MODE::MP;
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << "Library construction mode (--lc) should be one of " << ART_LIB_CONST_MODE_SE << ", " << ART_LIB_CONST_MODE_PE << ", " << ART_LIB_CONST_MODE_MP << ".");
            exit(EXIT_FAILURE);
        }
        stream = vm_.count("stream");
        sep_flag = vm_.count("sep_flag");

        seq_file = vm_["seq_file"].as<string>();
        read_len = vm_["read_len"].as<int>();
        max_indel = vm_["max_indel"].as<int>();
        max_indel = vm_["max_indel"].as<int>();
        ins_rate_1 = vm_["ins_rate_1"].as<double>();
        ins_rate_2 = vm_["ins_rate_2"].as<double>();
        del_rate_1 = vm_["del_rate_1"].as<double>();
        del_rate_2 = vm_["del_rate_2"].as<double>();
        fcov = vm_["fcov"].as<string>();
        pe_frag_dist_mean = vm_["pe_frag_dist_mean"].as<int>();
        pe_frag_dist_std_dev = vm_["pe_frag_dist_std_dev"].as<double>();
        min_qual = vm_["min_qual"].as<int>();
        max_qual = vm_["max_qual"].as<int>();
        q_shift_1 = vm_["q_shift_1"].as<int>();
        q_shift_2 = vm_["q_shift_2"].as<int>();
        qual_file_1 = vm_["qual_file_1"].as<string>();
        qual_file_2 = vm_["qual_file_2"].as<string>();
        id = vm_["id"].as<string>();
        parallel = vm_["parallel"].as<int>();
}

po::options_description option_parser(){
    OutputDispatcherFactory out_dispatcher_factory_;
    out_dispatcher_factory_.add(std::make_shared<FastqReadOutputFactory>());
    out_dispatcher_factory_.add(std::make_shared<BamReadOutputFactory>());
    out_dispatcher_factory_.add(std::make_shared<HeadlessBamReadOutputFactory>());
    po::options_description general_opts("General Options");
    general_opts.add_options()("help", "print out usage information");
    general_opts.add_options()("version", "display version info");

    po::options_description required_opts("Required Options");

    required_opts.add_options()("mode", po::value<std::string>()->default_value("wgs"),
                                R"(simulation mode, should be "wgs", "trans" or "template")");
    required_opts.add_options()("lc", po::value<std::string>()->default_value("se"),
                                R"(library construction mode, should be "se", "pe" or "mp")");

    required_opts.add_options()("seq_file", po::value<std::string>(),
                                "the filename of input reference genome, reference "
                                "transcriptome, or templates");
    required_opts.add_options()("fcov", po::value<std::string>(),
                                "the fold of read coverage to be simulated or number of reads/read pairs "
                                "generated for each sequence for simulating cDNA reads, or a float for "
                                "simulating WGS reads.");

    po::options_description art_opts("ART-specific options");
    art_opts.add_options()(
        "id", po::value<std::string>()->default_value("ART"), "the prefix identification tag for read ID");

    art_opts.add_options()(
        "qual_file_1", po::value<std::string>()->default_value(""), "the first-read quality profile");
    art_opts.add_options()("qual_file_2", po::value<std::string>()->default_value(""),
                           "the second-read quality profile. For PE/MP only.");
    art_opts.add_options()(
        "ins_rate_1", po::value<double>()->default_value(DEFAULT_INS_RATE_1), "the first-read insertion rate");
    art_opts.add_options()(
        "ins_rate_2", po::value<double>()->default_value(DEFAULT_INS_RATE_2), "the second-read insertion rate");
    art_opts.add_options()(
        "del_rate_1", po::value<double>()->default_value(DEFAULT_DEL_RATE_1), "the second-read deletion rate");
    art_opts.add_options()(
        "del_rate_2", po::value<double>()->default_value(DEFAULT_DEL_RATE_2), "the second-read deletion rate");
    art_opts.add_options()("sep_flag",
                           "use separate quality profiles for different bases. Default "
                           "is to use same quality profile regardless its position");
    art_opts.add_options()("max_indel", po::value<int>()->default_value(DEFAULT_MAX_INDEL),
                           "the maximum total number of insertion and deletion per read");
    art_opts.add_options()("read_len", po::value<int>(), "read length to be simulated");
    art_opts.add_options()("pe_frag_dist_mean", po::value<int>()->default_value(0),
                           "Mean distance between DNA/RNA fragments for paired-end simulations");
    art_opts.add_options()("pe_frag_dist_std_dev", po::value<double>()->default_value(0),
                           "Std. deviation of distance between DNA/RNA fragments for paired-end "
                           "simulations");
    art_opts.add_options()(
        "q_shift_1", po::value<int>()->default_value(0), "the amount to shift every first-read quality score by");
    art_opts.add_options()(
        "q_shift_2", po::value<int>()->default_value(0), "the amount to shift every second-read quality score by");
    art_opts.add_options()("min_qual", po::value<int>()->default_value(MIN_QUAL), "the minimum base quality score");
    art_opts.add_options()("max_qual", po::value<int>()->default_value(MAX_QUAL), "the maximum base quality score");

    po::options_description parallel_opts("Parallelism-related options");
    parallel_opts.add_options()("parallel", po::value<int>()->default_value(PARALLEL_ALL),
                                "Parallel level. -1 for disable, 0 for all CPUs, >=1 to specify number "
                                "of threads.");
    art_opts.add_options()("stream",
                           "If specified, will use streamline FASTA parser (For "
                           "transcriptome/amplicon); Otherwise will use HTSLib indexed "
                           "FASTA parser (For genome).");
    po::options_description po_desc;
    po_desc.add(general_opts).add(required_opts);
    out_dispatcher_factory_.patch_options(po_desc);
    po_desc.add(art_opts).add(parallel_opts);
    return po_desc;
}
} // art_modern
} // labw