#include "art/ArtCmdOpts.hh"

#include "art/ArtConstants.hh"
#include "art/ArtParams.hh"
#include "art/Empdist.hh"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/program_options.hpp>

#include <htslib/faidx.h>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

#include "CExceptionsProxy.hh"
#include "art_modern_constants.hh"

#include "ds/CoverageInfo.hh"

#include "out/OutputDispatcher.hh"

#include "utils/fs_utils.hh"
#include "utils/mpi_utils.hh"
#include "utils/param_utils.hh"
#include "utils/version_utils.hh"

namespace po = boost::program_options;

namespace labw::art_modern {
namespace {

constexpr char ARG_VERSION[] = "version";
constexpr char ARG_HELP[] = "help";
constexpr char ARG_SIMULATION_MODE[] = "mode";
constexpr char ARG_LIB_CONST_MODE[] = "lc";

constexpr char ARG_INPUT_FILE_NAME[] = "i-file";
constexpr char ARG_INPUT_FILE_PARSER[] = "i-parser";
constexpr char ARG_INPUT_FILE_TYPE[] = "i-type";
constexpr char ARG_FCOV[] = "i-fcov";
constexpr char ARG_BATCH_SIZE[] = "i-batch_size";

constexpr char ARG_ID[] = "id";
constexpr char ARG_PARALLEL[] = "parallel";
constexpr char ARG_QUAL_FILE_1[] = "qual_file_1";
constexpr char ARG_QUAL_FILE_2[] = "qual_file_2";
constexpr char ARG_READ_LEN[] = "read_len";
constexpr char ARG_MAX_INDEL[] = "max_indel";
constexpr char ARG_MAX_N[] = "max_n";
constexpr char ARG_INS_RATE_1[] = "ins_rate_1";
constexpr char ARG_INS_RATE_2[] = "ins_rate_2";
constexpr char ARG_DEL_RATE_1[] = "del_rate_1";
constexpr char ARG_DEL_RATE_2[] = "del_rate_2";
constexpr char ARG_SEP_FLAG[] = "sep_flag";
constexpr char ARG_PE_FRAG_DIST_MEAN[] = "pe_frag_dist_mean";
constexpr char ARG_PE_FRAG_DIST_STD_DEV[] = "pe_frag_dist_std_dev";
constexpr char ARG_MIN_QUAL[] = "min_qual";
constexpr char ARG_MAX_QUAL[] = "max_qual";
constexpr char ARG_Q_SHIFT_1[] = "q_shift_1";
constexpr char ARG_Q_SHIFT_2[] = "q_shift_2";

po::options_description option_parser() noexcept
{
    const OutputDispatcherFactory& out_dispatcher_factory_ = get_output_dispatcher_factory();
    po::options_description general_opts("General Options");
    general_opts.add_options()(static_cast<const char*>(ARG_HELP), "print out usage information");
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
            + INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS + ".")
            .c_str());
    required_opts.add_options()(ARG_BATCH_SIZE, po::value<int>()->default_value(DEFAULT_BATCH_SIZE),
        (std::string() + "Batch size for " + INPUT_FILE_PARSER_STREAM + " input parser").c_str());

    required_opts.add_options()(ARG_INPUT_FILE_NAME, po::value<std::string>(),
        "the filename of input reference genome, reference "
        "transcriptome, or templates");
    required_opts.add_options()(ARG_FCOV, po::value<std::string>()->default_value("0.0"),
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
    art_opts.add_options()(ARG_MAX_N, po::value<int>()->default_value(DEFAULT_MAX_N),
        "the maximum total number of ambiguous bases (N) per read");
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

    try {
        store(parse_command_line(argc, argv, po_desc), vm_);
        notify(vm_);
    } catch (const std::exception& exp) {
        BOOST_LOG_TRIVIAL(fatal) << exp.what();
        print_help(po_desc);
        abort_mpi();
    }

    if (vm_.count(ARG_VERSION) != 0U) {
        print_version();
        bye_mpi();
        exit_mpi(EXIT_SUCCESS);
    }
    if (vm_.count(ARG_HELP) != 0U) {
        print_help(po_desc);
        bye_mpi();
        exit_mpi(EXIT_SUCCESS);
    }
    return vm_;
}

SIMULATION_MODE get_simulation_mode(const std::string& simulation_mode_str)
{
    if (simulation_mode_str == SIMULATION_MODE_WGS) {
        return SIMULATION_MODE::WGS;
    }
    if (simulation_mode_str == SIMULATION_MODE_TRANS) {
        return SIMULATION_MODE::TRANS;
    }
    if (simulation_mode_str == SIMULATION_MODE_TEMPLATE) {
        return SIMULATION_MODE::TEMPLATE;
    }

    BOOST_LOG_TRIVIAL(fatal) << "Simulation mode (--" << ARG_SIMULATION_MODE << ") should be one of "
                             << SIMULATION_MODE_WGS << ", " << SIMULATION_MODE_TRANS << ", " << SIMULATION_MODE_TEMPLATE
                             << ".";
    abort_mpi();
}

ART_LIB_CONST_MODE get_art_lib_const_mode(const std::string& lib_const_mode_str)
{
    for (int i =0; i < 3;i++) {
        if (lib_const_mode_str == ART_LIB_CONST_MODE_STR[i]) {
            return static_cast<ART_LIB_CONST_MODE>(i);
        }
    }

    BOOST_LOG_TRIVIAL(fatal) << "Library construction mode (--" << ARG_LIB_CONST_MODE << ") should be one of "
                             << ART_LIB_CONST_MODE_SE << ", " << ART_LIB_CONST_MODE_PE << ", " << ART_LIB_CONST_MODE_MP
                             << ".";
    abort_mpi();
}

INPUT_FILE_TYPE get_input_file_type(const std::string& input_file_type_str, const std::string& input_file_name)
{
    if (input_file_type_str == INPUT_FILE_TYPE_FASTA) {
        return INPUT_FILE_TYPE::FASTA;
    }
    if (input_file_type_str == INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS) {
        return INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS;
    }
    if (input_file_type_str == INPUT_FILE_TYPE_AUTO) {
        for (const auto& fasta_file_end : std::vector<std::string> { ".fna", ".fsa", ".fa", ".fasta" }) {
            if (boost::algorithm::ends_with(input_file_name, fasta_file_end)) {
                return INPUT_FILE_TYPE::FASTA;
            }
        }
        BOOST_LOG_TRIVIAL(fatal) << "Automatic inference of input file type failed! Modify value of this param (--"
                                 << ARG_INPUT_FILE_TYPE << ") to be one of " << INPUT_FILE_TYPE_FASTA << ", "
                                 << INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS << ".";
        abort_mpi();
    }

    BOOST_LOG_TRIVIAL(fatal) << "Input file type (--" << ARG_INPUT_FILE_TYPE << ") should be one of "
                             << INPUT_FILE_TYPE_FASTA << ", " << INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS << ", "
                             << INPUT_FILE_TYPE_AUTO << ".";
    abort_mpi();
}

INPUT_FILE_PARSER get_input_file_parser(
    const std::string& input_file_parser_str, const std::string& input_file_path, const SIMULATION_MODE simulation_mode)
{
    if (input_file_parser_str == INPUT_FILE_PARSER_MEMORY) {
        return INPUT_FILE_PARSER::MEMORY;
    }
    if (input_file_parser_str == INPUT_FILE_PARSER_HTSLIB) {
        return INPUT_FILE_PARSER::HTSLIB;
    }
    if (input_file_parser_str == INPUT_FILE_PARSER_STREAM) {
        return INPUT_FILE_PARSER::STREAM;
    }
    if (input_file_parser_str == INPUT_FILE_PARSER_AUTO) {
        const auto file_size = get_file_size(input_file_path);
        const auto file_too_large = file_size == -1 || file_size > G_SIZE;
        if (simulation_mode == SIMULATION_MODE::WGS) {
            if (file_too_large) {
                return INPUT_FILE_PARSER::HTSLIB;
            }
            return INPUT_FILE_PARSER::MEMORY;
        }

        if (file_too_large) {
            return INPUT_FILE_PARSER::STREAM;
        }
        return INPUT_FILE_PARSER::MEMORY;
    }

    BOOST_LOG_TRIVIAL(fatal) << "Input file parser (--" << ARG_INPUT_FILE_PARSER << ") should be one of "
                             << INPUT_FILE_PARSER_MEMORY << ", " << INPUT_FILE_PARSER_HTSLIB << ", "
                             << INPUT_FILE_PARSER_STREAM << ", " << INPUT_FILE_PARSER_AUTO << ".";
    abort_mpi();
}

CoverageInfo get_coverage_info(
    const std::string& fcov_arg_str, const INPUT_FILE_TYPE input_file_type, const SIMULATION_MODE simulation_mode)
{
    if (input_file_type == INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS) {
        return CoverageInfo(0.0);
    }

    try {
        auto d = boost::lexical_cast<double>(fcov_arg_str);
        if (simulation_mode == SIMULATION_MODE::TEMPLATE) {
            return CoverageInfo(d, 0.0);
        }

        return CoverageInfo(d);

    } catch (const boost::bad_lexical_cast&) {
        validate_input_filename(fcov_arg_str, ARG_FCOV);
        std::ifstream cov_fs(fcov_arg_str, std::ios::binary);
        auto coverage_info = CoverageInfo(cov_fs);
        cov_fs.close();
        return coverage_info;
    }
}

void validate_min_max_qual(const int min_qual, const int max_qual)
{
    if (min_qual < 0 || min_qual > MAX_QUAL) {
        BOOST_LOG_TRIVIAL(fatal) << "Input Error: The minimum quality score must be an integer in [0," << MAX_QUAL
                                 << "]";
        abort_mpi();
    }
    if (max_qual <= min_qual || max_qual > MAX_QUAL) {
        BOOST_LOG_TRIVIAL(fatal) << "Input Error: The quality score must be an integer in [" << min_qual << ", "
                                 << MAX_QUAL << "]";
        abort_mpi();
    }
}

void validate_qual_files(
    const std::string& qual_file_1, const std::string& qual_file_2, const ART_LIB_CONST_MODE art_lib_const_mode)
{
    validate_input_filename(qual_file_1, ARG_QUAL_FILE_1);
    if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        validate_input_filename(qual_file_2, ARG_QUAL_FILE_2);
    }
}

Empdist read_emp(const std::string& qual_file_1, const std::string& qual_file_2, const size_t read_len,
    const ART_LIB_CONST_MODE art_lib_const_mode, const bool sep_flag, const int q_shift_1, const int q_shift_2,
    const int min_qual, const int max_qual)
{
    validate_min_max_qual(min_qual, max_qual);
    validate_qual_files(qual_file_1, qual_file_2, art_lib_const_mode);
    auto qdist = Empdist(qual_file_1, qual_file_2, sep_flag, art_lib_const_mode != ART_LIB_CONST_MODE::SE, read_len);
    size_t r1_profile_size = 0;
    size_t r2_profile_size = 0;
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
        abort_mpi();
    }

    if (read_len > r2_profile_size && art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        if (r2_profile_size == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: " << qual_file_2 << ", is not a valid profile.";
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length, " << read_len
                                     << ", exceeds the maximum second read profile length, " << r2_profile_size << ".";
        }
        abort_mpi();
    }
    qdist.shift_all_emp(sep_flag, q_shift_1, q_shift_2, min_qual, max_qual);
    return qdist;
}

void validate_htslib_parser(const std::string& input_file_path)
{
    const char* fasta_path = input_file_path.c_str();
    BOOST_LOG_TRIVIAL(info) << "HTSLib parser requested. Checking FAI...";
    const auto seq_file_fai_path
        = std::string(CExceptionsProxy::assert_not_null(fai_path(fasta_path), USED_HTSLIB_NAME, "Failed to load FAI"));
    if (!exists(boost::filesystem::path(seq_file_fai_path))) {
        BOOST_LOG_TRIVIAL(info) << "Building missing FAI...";
        CExceptionsProxy::assert_numeric(fai_build(fasta_path), USED_HTSLIB_NAME, "Failed to build FAI");
    } else {
        BOOST_LOG_TRIVIAL(info) << "Loading existing FAI...";
        CExceptionsProxy::assert_not_null(
            fai_load_format(fasta_path, FAI_FASTA), USED_HTSLIB_NAME, "Failed to load FAI");
    }
}

std::vector<double> gen_per_base_mutation_rate(const int read_len, const double p, const int max_indel)
{
    std::vector<double> rate;
    if (max_indel == 0 || p < 1E-30) {
        return rate;
    }

    double tp = 0;
    double p_cdf = 0;
    for (auto i = 0; i < read_len; i++) {
        tp = cdf(complement(boost::math::binomial(read_len, p), i));
        rate.emplace_back(tp);
        if (max_indel > 0 && i >= max_indel) {
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
        abort_mpi();
    }
    if (read_len == 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Read length must be specified.";
        abort_mpi();
    }
}

int validate_parallel(const int parallel_arg)
{
    int parallel = parallel_arg;
    const auto max_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (parallel_arg == PARALLEL_ALL) {
        parallel = max_threads;
    } else if (parallel_arg == PARALLEL_DISABLE) {
        parallel = 1;
    } else if (parallel_arg > max_threads) {
        BOOST_LOG_TRIVIAL(warning) << "parallel (" << parallel
                                   << ") is greater than the "
                                      "maximum number of threads available on the system ("
                                   << max_threads << ").";
    } else if (parallel_arg < -1) {
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
            if (pe_frag_dist_std_dev <= 0 || pe_frag_dist_mean <= 0) {
                BOOST_LOG_TRIVIAL(fatal) << "set pe_frag_dist_std_dev and "
                                            "pe_frag_dist_mean for PE reads for "
                                         << SIMULATION_MODE_WGS << " or " << SIMULATION_MODE_TRANS << " mode)";
                abort_mpi();
            }

            if (pe_frag_dist_mean <= read_len) {
                BOOST_LOG_TRIVIAL(fatal) << "The read length must be shorter than the "
                                            "pe_frag_dist_mean fragment length specified.";
                abort_mpi();
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
    if (input_file_type == INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS) {
        if (art_simulation_mode == SIMULATION_MODE::WGS) {
            BOOST_LOG_TRIVIAL(fatal) << "Input using PBSIM3 transcripts format are not supported for WGS simulation.";
            abort_mpi();
        }
        if (input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
            BOOST_LOG_TRIVIAL(fatal) << "Input using PBSIM3 transcripts format are not supported for HTSLib parser";
            abort_mpi();
        }
    }
    if (art_simulation_mode == SIMULATION_MODE::WGS && input_file_parser == INPUT_FILE_PARSER::STREAM) {
        BOOST_LOG_TRIVIAL(fatal) << "STREAM parser is not supported for WGS simulation ";
        abort_mpi();
    }
    if (art_simulation_mode != SIMULATION_MODE::WGS && input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
        BOOST_LOG_TRIVIAL(fatal) << "HTSLib parser is only supported in WGS simulation ";
        abort_mpi();
    }
}

} // namespace

ArtParams parse_args(const int argc, char** argv)
{
    const boost::program_options::options_description po_desc_ = option_parser();
    const OutputDispatcherFactory out_dispatcher_factory_ = get_output_dispatcher_factory();

    std::vector<std::string> args { argv, argv + argc };
    BOOST_LOG_TRIVIAL(info) << "ARGS: " << boost::algorithm::join(args, " ");

    const auto& vm_ = generate_vm_while_handling_help_version(po_desc_, argc, argv);
    const auto& art_simulation_mode = get_simulation_mode(get_param<std::string>(vm_, ARG_SIMULATION_MODE));
    const auto& art_lib_const_mode = get_art_lib_const_mode(get_param<std::string>(vm_, ARG_LIB_CONST_MODE));
    const auto& input_file_name = get_param<std::string>(vm_, ARG_INPUT_FILE_NAME);
    validate_input_filename(input_file_name, ARG_INPUT_FILE_NAME);
    const auto& input_file_type
        = get_input_file_type(get_param<std::string>(vm_, ARG_INPUT_FILE_TYPE), input_file_name);
    const auto& input_file_parser = get_input_file_parser(
        get_param<std::string>(vm_, ARG_INPUT_FILE_PARSER), input_file_name, art_simulation_mode);
    validate_comp_mtx(input_file_parser, art_simulation_mode, input_file_type);
    if (input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
        validate_htslib_parser(input_file_name);
    }
    auto coverage_info = get_coverage_info(get_param<std::string>(vm_, ARG_FCOV), input_file_type, art_simulation_mode);

    auto id = get_param<std::string>(vm_, ARG_ID);

    const auto sep_flag = vm_.count(ARG_SEP_FLAG) > 0;
    const auto max_indel = get_param<int>(vm_, ARG_MAX_INDEL);
    const auto max_n = get_param<int>(vm_, ARG_MAX_N);
    const auto read_len = get_param<int>(vm_, ARG_READ_LEN);
    validate_read_length(read_len);

    auto per_base_ins_rate_1 = gen_per_base_mutation_rate(read_len, get_param<double>(vm_, ARG_INS_RATE_1), max_indel);
    auto per_base_del_rate_1 = gen_per_base_mutation_rate(read_len, get_param<double>(vm_, ARG_DEL_RATE_1), max_indel);
    auto per_base_ins_rate_2 = gen_per_base_mutation_rate(read_len, get_param<double>(vm_, ARG_INS_RATE_2), max_indel);
    auto per_base_del_rate_2 = gen_per_base_mutation_rate(read_len, get_param<double>(vm_, ARG_DEL_RATE_2), max_indel);

    const auto pe_frag_dist_mean = get_param<double>(vm_, ARG_PE_FRAG_DIST_MEAN);
    const auto pe_frag_dist_std_dev = get_param<double>(vm_, ARG_PE_FRAG_DIST_STD_DEV);
    validate_pe_frag_dist(pe_frag_dist_mean, pe_frag_dist_std_dev, read_len, art_lib_const_mode, art_simulation_mode);
    const auto pe_dist_mean_minus_2_std = static_cast<hts_pos_t>(pe_frag_dist_mean - 2 * pe_frag_dist_std_dev);

    const auto& parallel = validate_parallel(get_param<int>(vm_, ARG_PARALLEL));

    auto qdist = read_emp(get_param<std::string>(vm_, ARG_QUAL_FILE_1), get_param<std::string>(vm_, ARG_QUAL_FILE_2),
        read_len, art_lib_const_mode, sep_flag, get_param<int>(vm_, ARG_Q_SHIFT_1), get_param<int>(vm_, ARG_Q_SHIFT_2),
        get_param<int>(vm_, ARG_MIN_QUAL), get_param<int>(vm_, ARG_MAX_QUAL));
    std::array<double, HIGHEST_QUAL> err_prob {};
    for (int i = 0; i < HIGHEST_QUAL; i++) {
        err_prob[i] = std::pow(10, -i / 10.0);
    }
    const auto& batch_size = get_param<int>(vm_, ARG_BATCH_SIZE);
    if (batch_size < 1) {
        BOOST_LOG_TRIVIAL(fatal) << "Batch size (" << batch_size << ") must be greater than 1";
        abort_mpi();
    }
    return { art_simulation_mode, art_lib_const_mode, input_file_name, input_file_type, input_file_parser, parallel,
        sep_flag, std::move(id), std::move(coverage_info), max_n, read_len, pe_frag_dist_mean, pe_frag_dist_std_dev,
        std::move(per_base_ins_rate_1), std::move(per_base_del_rate_1), std::move(per_base_ins_rate_2),
        std::move(per_base_del_rate_2), err_prob, pe_dist_mean_minus_2_std, std::move(qdist), batch_size, vm_,
        std::move(args) };
}

} // namespace labw::art_modern
