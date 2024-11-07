#include "ArtCmdOpts.hh"
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>

#include "ArtConstants.hh"
#include "CExceptionsProxy.hh"
#include "art_modern_constants.hh"
#include "fasta/FaidxFetch.hh"
#include "fasta/InMemoryFastaFetch.hh"
#include "fasta/Pbsim3TranscriptBatcher.hh"
#include "global_variables.hh"
#include "out/BamReadOutput.hh"
#include "out/FastqReadOutput.hh"
#include "out/HeadlessBamReadOutput.hh"
#include "out/OutputDispatcher.hh"
#include "utils/version_utils.hh"

namespace po = boost::program_options;

const char ARG_VERSION[] = "version";
const char ARG_HELP[] = "help";
const char ARG_SIMULATION_MODE[] = "mode";
const char ARG_LIB_CONST_MODE[] = "lc";

const char ARG_INPUT_FILE_NAME[] = "i-file";
const char ARG_INPUT_FILE_PARSER[] = "i-parser";
const char ARG_INPUT_FILE_TYPE[] = "i-type";
const char ARG_FCOV[] = "i-fcov";
const char ARG_BATCH_SIZE[] = "i-batch_size";

const char ARG_ID[] = "id";
const char ARG_PARALLEL[] = "parallel";
const char ARG_QUAL_FILE_1[] = "qual_file_1";
const char ARG_QUAL_FILE_2[] = "qual_file_2";
const char ARG_READ_LEN[] = "read_len";
const char ARG_MAX_INDEL[] = "max_indel";
const char ARG_INS_RATE_1[] = "ins_rate_1";
const char ARG_INS_RATE_2[] = "ins_rate_2";
const char ARG_DEL_RATE_1[] = "del_rate_1";
const char ARG_DEL_RATE_2[] = "del_rate_2";
const char ARG_SEP_FLAG[] = "sep_flag";
const char ARG_PE_FRAG_DIST_MEAN[] = "pe_frag_dist_mean";
const char ARG_PE_FRAG_DIST_STD_DEV[] = "pe_frag_dist_std_dev";
const char ARG_MIN_QUAL[] = "min_qual";
const char ARG_MAX_QUAL[] = "max_qual";
const char ARG_Q_SHIFT_1[] = "q_shift_1";
const char ARG_Q_SHIFT_2[] = "q_shift_2";

namespace labw::art_modern {

OutputDispatcherFactory get_output_dispatcher_factory() noexcept
{
    OutputDispatcherFactory out_dispatcher_factory;
    out_dispatcher_factory.add(new FastqReadOutputFactory());
    out_dispatcher_factory.add(new BamReadOutputFactory());
    out_dispatcher_factory.add(new HeadlessBamReadOutputFactory());
    return out_dispatcher_factory;
}

po::options_description option_parser() noexcept
{
    OutputDispatcherFactory out_dispatcher_factory_ = get_output_dispatcher_factory();
    po::options_description general_opts("General Options");
    general_opts.add_options()(ARG_HELP, "print out usage information");
    general_opts.add_options()(ARG_VERSION, "display version info");

    po::options_description required_opts("Required Options");

    required_opts.add_options()(ARG_SIMULATION_MODE, po::value<std::string>()->default_value(SIMULATION_MODE_WGS),
        (std::string() + "simulation mode, should be " + SIMULATION_MODE_WGS + ", " + SIMULATION_MODE_TRANS + ", "
            + SIMULATION_MODE_TEMPLATE + ".")
            .c_str());
    required_opts.add_options()(ARG_LIB_CONST_MODE, po::value<std::string>()->default_value(ART_LIB_CONST_MODE_SE),
        (std::string() + "library construction mode, should be " + ART_LIB_CONST_MODE_SE + ", " + ART_LIB_CONST_MODE_PE
            + ", " + ART_LIB_CONST_MODE_MP + ".")
            .c_str());
    required_opts.add_options()(ARG_INPUT_FILE_PARSER, po::value<std::string>()->default_value(INPUT_FILE_PARSER_AUTO),
        (std::string() + "input file parser, should be " + INPUT_FILE_PARSER_AUTO + ", " + INPUT_FILE_PARSER_MEMORY
            + ", " + INPUT_FILE_PARSER_HTSLIB + ", " + INPUT_FILE_PARSER_STREAM + ".")
            .c_str());
    required_opts.add_options()(ARG_INPUT_FILE_TYPE, po::value<std::string>()->default_value(INPUT_FILE_TYPE_AUTO),
        (std::string() + "input file type, should be " + INPUT_FILE_TYPE_AUTO + ", " + INPUT_FILE_TYPE_FASTA + ", "
            + INPUT_FILE_TYPE_PBSIM3_TEMPLATE + ".")
            .c_str());
    required_opts.add_options()(ARG_BATCH_SIZE, po::value<int>()->default_value(DEFAULT_BATCH_SIZE),
        (std::string() + "Batch size for " + INPUT_FILE_PARSER_STREAM + " input parser").c_str());

    required_opts.add_options()(ARG_INPUT_FILE_NAME, po::value<std::string>(),
        "the filename of input reference genome, reference "
        "transcriptome, or templates");
    required_opts.add_options()(ARG_FCOV, po::value<std::string>(),
        "the fold of read coverage to be simulated or number of reads/read pairs "
        "generated for each sequence for simulating cDNA reads, or a double for "
        "simulating WGS reads.");

    po::options_description art_opts("ART-specific options");
    art_opts.add_options()(
        ARG_ID, po::value<std::string>()->default_value(ART_PROGRAM_NAME), "the prefix identification tag for read ID");

    art_opts.add_options()(
        ARG_QUAL_FILE_1, po::value<std::string>()->default_value(""), "the first-read quality profile");
    art_opts.add_options()(ARG_QUAL_FILE_2, po::value<std::string>()->default_value(""),
        "the second-read quality profile. For PE/MP only.");
    art_opts.add_options()(
        ARG_INS_RATE_1, po::value<double>()->default_value(DEFAULT_INS_RATE_1), "the first-read insertion rate");
    art_opts.add_options()(
        ARG_INS_RATE_2, po::value<double>()->default_value(DEFAULT_INS_RATE_2), "the second-read insertion rate");
    art_opts.add_options()(
        ARG_DEL_RATE_1, po::value<double>()->default_value(DEFAULT_DEL_RATE_1), "the second-read deletion rate");
    art_opts.add_options()(
        ARG_DEL_RATE_2, po::value<double>()->default_value(DEFAULT_DEL_RATE_2), "the second-read deletion rate");
    art_opts.add_options()(ARG_SEP_FLAG,
        "use separate quality profiles for different bases. Default "
        "is to use same quality profile regardless its position");
    art_opts.add_options()(ARG_MAX_INDEL, po::value<int>()->default_value(DEFAULT_MAX_INDEL),
        "the maximum total number of insertion and deletion per read");
    art_opts.add_options()(ARG_READ_LEN, po::value<int>(), "read length to be simulated");
    art_opts.add_options()(ARG_PE_FRAG_DIST_MEAN, po::value<double>()->default_value(0),
        "Mean distance between DNA/RNA fragments for paired-end simulations");
    art_opts.add_options()(ARG_PE_FRAG_DIST_STD_DEV, po::value<double>()->default_value(0),
        "Std. deviation of distance between DNA/RNA fragments for paired-end "
        "simulations");
    art_opts.add_options()(
        ARG_Q_SHIFT_1, po::value<int>()->default_value(0), "the amount to shift every first-read quality score by");
    art_opts.add_options()(
        ARG_Q_SHIFT_2, po::value<int>()->default_value(0), "the amount to shift every second-read quality score by");
    art_opts.add_options()(ARG_MIN_QUAL, po::value<int>()->default_value(MIN_QUAL), "the minimum base quality score");
    art_opts.add_options()(ARG_MAX_QUAL, po::value<int>()->default_value(MAX_QUAL), "the maximum base quality score");

    po::options_description parallel_opts("Parallelism-related options");
    parallel_opts.add_options()(ARG_PARALLEL, po::value<int>()->default_value(PARALLEL_ALL),
        "Parallel level. -1 for disable, 0 for all CPUs, >=1 to specify number "
        "of threads.");
    po::options_description po_desc;
    po_desc.add(general_opts).add(required_opts);
    out_dispatcher_factory_.patch_options(po_desc);
    po_desc.add(art_opts).add(parallel_opts);
    return po_desc;
}

void print_help(const po::options_description& po_desc) { std::cout << po_desc << std::endl; }

po::variables_map generate_vm_while_handling_help_version(
    const po::options_description& po_desc, const int argc, char** argv)
{
    po::variables_map vm_;
    BOOST_LOG_TRIVIAL(info) << "ARGS: " << boost::algorithm::join(args, " ");

    try {
        store(po::parse_command_line(argc, argv, po_desc), vm_);
        notify(vm_);
    } catch (const std::exception& exp) {
        BOOST_LOG_TRIVIAL(fatal) << exp.what();
        print_help(po_desc);
        throw ArtCmdException();
    }

    if (vm_.count(ARG_VERSION)) {
        print_version();
        throw ArtCmdNormalExit();
    } else if (vm_.count(ARG_HELP)) {
        print_help(po_desc);
        throw ArtCmdNormalExit();
    }
    return vm_;
}

SIMULATION_MODE get_simulation_mode(const std::string& simulation_mode_str)
{
    if (simulation_mode_str == SIMULATION_MODE_WGS) {
        return SIMULATION_MODE::WGS;
    } else if (simulation_mode_str == SIMULATION_MODE_TRANS) {
        return SIMULATION_MODE::TRANS;
    } else if (simulation_mode_str == SIMULATION_MODE_TEMPLATE) {
        return SIMULATION_MODE::TEMPLATE;
    } else {
        BOOST_LOG_TRIVIAL(fatal) << "Simulation mode (--" << ARG_SIMULATION_MODE << ") should be one of "
                                 << SIMULATION_MODE_WGS << ", " << SIMULATION_MODE_TRANS << ", "
                                 << SIMULATION_MODE_TEMPLATE << ".";
        throw ArtCmdException();
    }
}

ART_LIB_CONST_MODE get_art_lib_const_mode(const std::string& lib_const_mode_str)
{
    if (lib_const_mode_str == ART_LIB_CONST_MODE_SE) {
        return ART_LIB_CONST_MODE::SE;
    } else if (lib_const_mode_str == ART_LIB_CONST_MODE_PE) {
        return ART_LIB_CONST_MODE::PE;
    } else if (lib_const_mode_str == ART_LIB_CONST_MODE_MP) {
        return ART_LIB_CONST_MODE::MP;
    } else {
        BOOST_LOG_TRIVIAL(fatal) << "Library construction mode (--" << ARG_LIB_CONST_MODE << ") should be one of "
                                 << ART_LIB_CONST_MODE_SE << ", " << ART_LIB_CONST_MODE_PE << ", "
                                 << ART_LIB_CONST_MODE_MP << ".";
        throw ArtCmdException();
    }
}

INPUT_FILE_TYPE get_input_file_type(const std::string& input_file_type_str, const std::string& input_file_name)
{
    if (input_file_type_str == INPUT_FILE_TYPE_FASTA) {
        return INPUT_FILE_TYPE::FASTA;
    } else if (input_file_type_str == INPUT_FILE_TYPE_PBSIM3_TEMPLATE) {
        return INPUT_FILE_TYPE::PBSIM3_TEMPLATE;
    } else if (input_file_type_str == INPUT_FILE_TYPE_AUTO) {
        for (const auto& fasta_file_end : std::vector<std::string> { "fna", "faa", "fa", "fasta" }) {
            if (boost::algorithm::ends_with(input_file_name, fasta_file_end)) {
                return INPUT_FILE_TYPE::FASTA;
            }
        }
        BOOST_LOG_TRIVIAL(fatal) << "Automatic inference of input file type failed! Modify value of this param (--"
                                 << ARG_INPUT_FILE_TYPE << ") to be one of " << INPUT_FILE_TYPE_FASTA << ", "
                                 << INPUT_FILE_TYPE_PBSIM3_TEMPLATE << ".";
        throw ArtCmdException();
    } else {
        BOOST_LOG_TRIVIAL(fatal) << "Input file type (--" << ARG_INPUT_FILE_TYPE << ") should be one of "
                                 << INPUT_FILE_TYPE_FASTA << ", " << INPUT_FILE_TYPE_PBSIM3_TEMPLATE << ", "
                                 << INPUT_FILE_TYPE_AUTO << ".";
        throw ArtCmdException();
    }
}

long get_file_size(const std::string& file_path) noexcept
{
    if (!boost::filesystem::is_regular_file(file_path)) {
        return -1;
    } else {
        try {
            return static_cast<long>(boost::filesystem::file_size(file_path));
        } catch (const boost::filesystem::filesystem_error&) {
            return -1;
        }
    }
}

INPUT_FILE_PARSER get_input_file_parser(
    const std::string& input_file_parser_str, const std::string& input_file_path, const SIMULATION_MODE simulation_mode)
{
    if (input_file_parser_str == INPUT_FILE_PARSER_MEMORY) {
        return INPUT_FILE_PARSER::MEMORY;
    } else if (input_file_parser_str == INPUT_FILE_PARSER_HTSLIB) {
        return INPUT_FILE_PARSER::HTSLIB;
    } else if (input_file_parser_str == INPUT_FILE_PARSER_STREAM) {
        return INPUT_FILE_PARSER::STREAM;
    } else if (input_file_parser_str == INPUT_FILE_PARSER_AUTO) {
        auto file_size = get_file_size(input_file_path);
        auto file_too_large = file_size == -1 || file_size > (1 * 1024 * 1024 * 1024);
        if (simulation_mode == SIMULATION_MODE::WGS) {
            if (file_too_large) {
                return INPUT_FILE_PARSER::HTSLIB;
            }
            return INPUT_FILE_PARSER::MEMORY;
        } else {
            if (file_too_large) {
                return INPUT_FILE_PARSER::STREAM;
            }
            return INPUT_FILE_PARSER::MEMORY;
        }
    } else {
        BOOST_LOG_TRIVIAL(fatal) << "Input file parser (--" << ARG_INPUT_FILE_PARSER << ") should be one of "
                                 << INPUT_FILE_PARSER_MEMORY << ", " << INPUT_FILE_PARSER_HTSLIB << ", "
                                 << INPUT_FILE_PARSER_STREAM << ", " << INPUT_FILE_PARSER_AUTO << ".";
        throw ArtCmdException();
    }
}

std::pair<CoverageInfo, BaseFastaFetch*> get_coverage_info_fasta_fetch(const std::string& fcov_arg_str,
    const INPUT_FILE_TYPE input_file_type, const INPUT_FILE_PARSER input_file_parser,
    const SIMULATION_MODE simulation_mode, const std::string& input_file_name)
{
    BaseFastaFetch* fasta_fetch;
    if (input_file_type == INPUT_FILE_TYPE::PBSIM3_TEMPLATE) {
        if (input_file_parser == INPUT_FILE_PARSER::MEMORY) {
            std::ifstream input_file_stream(input_file_name);
            Pbsim3TranscriptBatcher batcher(std::numeric_limits<int>::max(), input_file_stream);
            fasta_fetch = new InMemoryFastaFetch(batcher.fetch().first);
            auto coverage_info = batcher.fetch().second;
            input_file_stream.close();
            return { coverage_info, fasta_fetch };
        } else if (input_file_parser == INPUT_FILE_PARSER::STREAM) {
            return { CoverageInfo(0.0), new InMemoryFastaFetch() };
        }
    }
    if (fcov_arg_str.empty()) {
        BOOST_LOG_TRIVIAL(fatal) << "Coverage parameter (--" << ARG_FCOV << ") is required.";
        throw ArtCmdException();
    }
    if (input_file_parser == INPUT_FILE_PARSER::STREAM) {
        fasta_fetch = new InMemoryFastaFetch(); // Empty
    } else if (input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
        fasta_fetch = new FaidxFetch(input_file_name);
    } else if (input_file_parser == INPUT_FILE_PARSER::MEMORY) {
        fasta_fetch = new InMemoryFastaFetch(input_file_name);
    }
    try {
        auto d = boost::lexical_cast<double>(fcov_arg_str);
        if (simulation_mode == SIMULATION_MODE::TEMPLATE) {
            auto coverage_info = CoverageInfo(d, 0.0);
            return { coverage_info, fasta_fetch };
        } else {
            auto coverage_info = CoverageInfo(d);
            return { coverage_info, fasta_fetch };
        }
    } catch (const boost::bad_lexical_cast&) {
        std::ifstream X_FOLD(fcov_arg_str, std::ios::binary);
        auto coverage_info = CoverageInfo(X_FOLD);
        X_FOLD.close();
        return { coverage_info, fasta_fetch };
    }
}

void shift_emp(std::map<int, int> map_to_process, const int q_shift, const int min_qual, const int max_qual)
{
    for (auto& map_to_proces : map_to_process) {
        if (q_shift != 0) {
            if (q_shift < 0 && (-q_shift > map_to_proces.second)) {
                map_to_proces.second = min_qual;
            } else {
                map_to_proces.second = std::min(map_to_proces.second + q_shift, max_qual);
            }
        }
        map_to_proces.second = std::min(std::max(map_to_proces.second, min_qual), max_qual);
    }
}

void validate_min_max_qual(const int min_qual, const int max_qual)
{
    if (min_qual < 0 || min_qual > MAX_QUAL) {
        BOOST_LOG_TRIVIAL(fatal) << "Input Error: The minimum quality score must be an integer in [0," << MAX_QUAL
                                 << "]";
        throw ArtCmdException();
    }
    if (max_qual <= min_qual || max_qual > MAX_QUAL) {
        BOOST_LOG_TRIVIAL(fatal) << "Input Error: The quality score must be an integer in [" << min_qual << ", "
                                 << MAX_QUAL << "]";
        throw ArtCmdException();
    }
}

void shift_all_emp(const Empdist& qdist, const bool sep_flag, const int q_shift_1, const int q_shift_2,
    const int min_qual, const int max_qual)
{
    if (!sep_flag) {
        for (const auto& i : qdist.qual_dist_first) {
            shift_emp(i, q_shift_1, min_qual, max_qual);
        }
        for (const auto& i : qdist.qual_dist_second) {
            shift_emp(i, q_shift_2, min_qual, max_qual);
        }
    } else {
        for (const auto& i : qdist.a_qual_dist_first) {
            shift_emp(i, q_shift_1, min_qual, max_qual);
        }
        for (const auto& i : qdist.a_qual_dist_second) {
            shift_emp(i, q_shift_2, min_qual, max_qual);
        }

        for (const auto& i : qdist.c_qual_dist_first) {
            shift_emp(i, q_shift_1, min_qual, max_qual);
        }
        for (const auto& i : qdist.c_qual_dist_second) {
            shift_emp(i, q_shift_2, min_qual, max_qual);
        }

        for (const auto& i : qdist.g_qual_dist_first) {
            shift_emp(i, q_shift_1, min_qual, max_qual);
        }
        for (const auto& i : qdist.g_qual_dist_second) {
            shift_emp(i, q_shift_2, min_qual, max_qual);
        }

        for (const auto& i : qdist.t_qual_dist_first) {
            shift_emp(i, q_shift_1, min_qual, max_qual);
        }
        for (const auto& i : qdist.t_qual_dist_second) {
            shift_emp(i, q_shift_2, min_qual, max_qual);
        }
    }
}

void validate_input_filename(const std::string& input_file_path, const std::string& arg_name)
{
    if (input_file_path.empty()) {
        BOOST_LOG_TRIVIAL(fatal) << "An input file path for --" << arg_name << " must be specified.";
        throw ArtCmdException();
    }
    if (!boost::filesystem::exists(input_file_path)) {
        BOOST_LOG_TRIVIAL(fatal) << "Input file for --" << arg_name << " at '" << input_file_path
                                 << "' does not exist.";
        throw ArtCmdException();
    }
    if (!boost::filesystem::is_regular_file(input_file_path)) {
        BOOST_LOG_TRIVIAL(warning) << "Input file for --" << arg_name << " at '" << input_file_path
                                   << "' is not a regular file.";
    }
}
void validate_qual_files(
    const std::string& qual_file_1, const std::string& qual_file_2, const ART_LIB_CONST_MODE art_lib_const_mode)
{
    validate_input_filename(qual_file_1, ARG_QUAL_FILE_1);
    if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        validate_input_filename(qual_file_2, ARG_QUAL_FILE_1);
    }
}

Empdist read_emp(const std::string& qual_file_1, const std::string& qual_file_2, const size_t read_len,
    const ART_LIB_CONST_MODE art_lib_const_mode, const bool sep_flag, const int q_shift_1, const int q_shift_2,
    const int min_qual, const int max_qual)
{
    validate_min_max_qual(min_qual, max_qual);
    validate_qual_files(qual_file_1, qual_file_2, art_lib_const_mode);
    auto qdist = Empdist(qual_file_1, qual_file_2, sep_flag);
    size_t r1_profile_size;
    size_t r2_profile_size;
    if (sep_flag) {
        r1_profile_size = qdist.a_qual_dist_first.size();
        r2_profile_size = qdist.a_qual_dist_second.size();
    } else {
        r1_profile_size = qdist.qual_dist_first.size();
        r2_profile_size = qdist.qual_dist_second.size();
    }

    if (read_len > r1_profile_size) {
        if (r1_profile_size == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: " << qual_file_1 << ", is not a valid profile.";
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length, " << read_len
                                     << ", exceeds the maximum first read profile length, " << r1_profile_size << ".";
        }
        throw ArtCmdException();
    }

    if ((read_len > r2_profile_size) && art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        if (r2_profile_size == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: " << qual_file_2 << ", is not a valid profile.";
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length, " << read_len
                                     << ", exceeds the maximum second read profile length, " << r2_profile_size << ".";
        }
        throw ArtCmdException();
    }
    if (q_shift_1 != 0 || q_shift_2 != 0) {
        shift_all_emp(qdist, sep_flag, q_shift_1, q_shift_2, min_qual, max_qual);
    }
    return qdist;
}

void validate_htslib_parser(const std::string& input_file_path)
{
    const char* fasta_path = input_file_path.c_str();
    BOOST_LOG_TRIVIAL(info) << "HTSLib parser requested. Checking FAI...";
    auto seq_file_fai_path = std::string(fai_path(fasta_path));
    if (!boost::filesystem::exists(boost::filesystem::path(seq_file_fai_path))) {
        BOOST_LOG_TRIVIAL(info) << "Building missing FAI...";
        CExceptionsProxy::assert_numeric(fai_build(fasta_path), USED_HTSLIB_NAME, "Failed to build FAI");
    } else {
        BOOST_LOG_TRIVIAL(info) << "Loading existing FAI...";
        CExceptionsProxy::assert_not_null(
            fai_load_format(fasta_path, FAI_FASTA), USED_HTSLIB_NAME, "Failed to load FAI");
    }
}

std::vector<double> gen_per_base_mutation_rate(const int read_len, const double p, const int max_num)
{
    std::vector<double> rate;
    if (max_num == 0 || p < 1E-30) {
        return rate;
    }

    double tp;
    double p_cdf = 0;
    for (auto i = 0; i < read_len; i++) {
        tp = boost::math::cdf(boost::math::complement(boost::math::binomial(read_len, p), i));
        rate.emplace_back(tp);
        if (max_num > 0 && (i >= max_num)) {
            break;
        }
        p_cdf += tp;
        if (p_cdf >= 0.999999) {
            break;
        }
    }
    return rate;
}

void validate_read_length(const int read_len)
{
    if (read_len < 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length must be a positive integer.";
        throw ArtCmdException();
    }
    if (read_len == 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Read length must be specified.";
        throw ArtCmdException();
    }
}

int validate_parallel(int parallel)
{
    auto max_threads = static_cast<int>(boost::thread::hardware_concurrency());
    if (parallel == PARALLEL_ALL) {
        parallel = max_threads;
    } else if (parallel == PARALLEL_DISABLE) {
        parallel = 1;
    } else if (parallel > max_threads) {
        BOOST_LOG_TRIVIAL(warning) << "parallel (" << parallel
                                   << ") is greater than the "
                                      "maximum number of threads available on the system ("
                                   << max_threads << ").";
    } else if (parallel < -1) {
        BOOST_LOG_TRIVIAL(fatal) << "parallel (" << parallel << ") must be greater than or equal to -1.";
    }
    return parallel;
}

void validate_pe_frag_dist(const double pe_frag_dist_mean, const double pe_frag_dist_std_dev, const int read_len,
    const ART_LIB_CONST_MODE art_lib_const_mode, const SIMULATION_MODE art_simulation_mode)
{
    if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        if (art_simulation_mode == SIMULATION_MODE::TEMPLATE) {
            if (pe_frag_dist_std_dev != 0 || pe_frag_dist_mean != 0) {
                BOOST_LOG_TRIVIAL(warning) << "pe_frag_dist_std_dev and "
                                              "pe_frag_dist_mean ignored for "
                                           << SIMULATION_MODE_TEMPLATE << " mode.";
            }
        } else {
            if (!(pe_frag_dist_std_dev > 0 && pe_frag_dist_mean > 0)) {
                BOOST_LOG_TRIVIAL(fatal) << "set pe_frag_dist_std_dev and "
                                            "pe_frag_dist_mean for PE reads for "
                                         << SIMULATION_MODE_WGS << " or " << SIMULATION_MODE_TRANS << " mode)";
                throw ArtCmdException();
            }

            if (pe_frag_dist_mean <= read_len) {
                BOOST_LOG_TRIVIAL(fatal) << "The read length must be shorter than the "
                                            "pe_frag_dist_mean fragment length specified.";
                throw ArtCmdException();
            }
        }
    } else {
        if (art_simulation_mode == SIMULATION_MODE::TEMPLATE && (pe_frag_dist_std_dev != 0 || pe_frag_dist_mean != 0)) {
            BOOST_LOG_TRIVIAL(warning) << "pe_frag_dist_std_dev and "
                                          "pe_frag_dist_mean ignored for "
                                       << ART_LIB_CONST_MODE_SE << " mode.";
        }
    }
}

void validate_comp_mtx(const INPUT_FILE_PARSER input_file_parser, const SIMULATION_MODE art_simulation_mode,
    const INPUT_FILE_TYPE input_file_type)
{
    if (input_file_type == INPUT_FILE_TYPE::PBSIM3_TEMPLATE) {
        if (art_simulation_mode == SIMULATION_MODE::WGS) {
            BOOST_LOG_TRIVIAL(fatal) << "Input using PBSim3 transcripts format are not supported for WGS simulation.";
            throw ArtCmdException();
        } else if (input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
            BOOST_LOG_TRIVIAL(fatal) << "Input using PBSim3 transcripts format are not supported for HTSLib parser";
            throw ArtCmdException();
        }
    }
    if (art_simulation_mode == SIMULATION_MODE::WGS && input_file_parser == INPUT_FILE_PARSER::STREAM) {
        BOOST_LOG_TRIVIAL(fatal) << "STREAM parser is not supported for WGS simulation ";
        throw ArtCmdException();
    }
    if (art_simulation_mode != SIMULATION_MODE::WGS && input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
        BOOST_LOG_TRIVIAL(fatal) << "HTSLib parser is only supported in WGS simulation ";
        throw ArtCmdException();
    }
}

ArtParams parse_args(int argc, char** argv)
{
    const boost::program_options::options_description po_desc_ = option_parser();
    const OutputDispatcherFactory out_dispatcher_factory_ = get_output_dispatcher_factory();

    for (int i = 0; i < argc; i++) {
        args.emplace_back(argv[i]);
    }
    auto vm_ = generate_vm_while_handling_help_version(po_desc_, argc, argv);
    auto art_simulation_mode = get_simulation_mode(vm_[ARG_SIMULATION_MODE].as<std::string>());
    auto art_lib_const_mode = get_art_lib_const_mode(vm_[ARG_LIB_CONST_MODE].as<std::string>());
    auto input_file_name = vm_[ARG_INPUT_FILE_NAME].as<std::string>();
    validate_input_filename(input_file_name, ARG_INPUT_FILE_NAME);
    auto input_file_type = get_input_file_type(vm_[ARG_INPUT_FILE_TYPE].as<std::string>(), input_file_name);
    auto input_file_parser
        = get_input_file_parser(vm_[ARG_INPUT_FILE_PARSER].as<std::string>(), input_file_name, art_simulation_mode);
    validate_comp_mtx(input_file_parser, art_simulation_mode, input_file_type);
    if (input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
        validate_htslib_parser(input_file_name);
    }
    auto ci_ff = get_coverage_info_fasta_fetch(
        vm_[ARG_FCOV].as<std::string>(), input_file_type, input_file_parser, art_simulation_mode, input_file_name);
    auto const& coverage_info = ci_ff.first;
    auto fasta_fetch = ci_ff.second;

    auto id = vm_[ARG_ID].as<std::string>();

    auto sep_flag = vm_.count(ARG_SEP_FLAG) > 0;
    auto max_indel = vm_[ARG_MAX_INDEL].as<int>();
    auto read_len = vm_[ARG_READ_LEN].as<int>();
    validate_read_length(read_len);

    auto per_base_ins_rate_1 = gen_per_base_mutation_rate(read_len, vm_[ARG_INS_RATE_1].as<double>(), max_indel);
    auto per_base_del_rate_1 = gen_per_base_mutation_rate(read_len, vm_[ARG_DEL_RATE_1].as<double>(), max_indel);
    auto per_base_ins_rate_2 = gen_per_base_mutation_rate(read_len, vm_[ARG_INS_RATE_2].as<double>(), max_indel);
    auto per_base_del_rate_2 = gen_per_base_mutation_rate(read_len, vm_[ARG_DEL_RATE_2].as<double>(), max_indel);

    auto pe_frag_dist_mean = vm_[ARG_PE_FRAG_DIST_MEAN].as<double>();
    auto pe_frag_dist_std_dev = vm_[ARG_PE_FRAG_DIST_STD_DEV].as<double>();
    validate_pe_frag_dist(pe_frag_dist_mean, pe_frag_dist_std_dev, read_len, art_lib_const_mode, art_simulation_mode);
    auto pe_dist_mean_minus_2_std = static_cast<hts_pos_t>(pe_frag_dist_mean - 2 * pe_frag_dist_std_dev);

    auto parallel = validate_parallel(vm_[ARG_PARALLEL].as<int>());

    auto qdist = read_emp(vm_[ARG_QUAL_FILE_1].as<std::string>(), vm_[ARG_QUAL_FILE_2].as<std::string>(), read_len,
        art_lib_const_mode, sep_flag, vm_[ARG_Q_SHIFT_1].as<int>(), vm_[ARG_Q_SHIFT_2].as<int>(),
        vm_[ARG_MIN_QUAL].as<int>(), vm_[ARG_MAX_QUAL].as<int>());
    std::array<double, HIGHEST_QUAL> err_prob {};

    for (int i = 0; i < HIGHEST_QUAL; i++) {
        err_prob[i] = std::pow(10, -i / 10.0);
    }
    auto batch_size = vm_[ARG_BATCH_SIZE].as<int>();
    if (batch_size < 1) {
        BOOST_LOG_TRIVIAL(fatal) << "Batch size (" << batch_size << ") must be greater than 1";
        throw ArtCmdException();
    }
    return { art_simulation_mode, art_lib_const_mode, input_file_name, input_file_type, input_file_parser, parallel,
        sep_flag, id, coverage_info, read_len, pe_frag_dist_mean, pe_frag_dist_std_dev, per_base_ins_rate_1,
        per_base_del_rate_1, per_base_ins_rate_2, per_base_del_rate_2, err_prob, pe_dist_mean_minus_2_std, qdist,
        batch_size, fasta_fetch, out_dispatcher_factory_.create(vm_, fasta_fetch) };
}

}