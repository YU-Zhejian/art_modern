#include "ArtCmdOpts.hh"
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>

#include "ArtConstants.hh"
#include "ArtVersion.hh"
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

using namespace std;
namespace po = boost::program_options;

#define ARG_VERSION "version"
#define ARG_HELP "help"
#define ARG_SIMULATION_MODE "mode"
#define ARG_LIB_CONST_MODE "lc"

#define ARG_INPUT_FILE_NAME "i-file"
#define ARG_INPUT_FILE_PARSER "i-parser"
#define ARG_INPUT_FILE_TYPE "i-type"
#define ARG_FCOV "i-fcov"

#define ARG_ID "id"
#define ARG_PARALLEL "parallel"
#define ARG_QUAL_FILE_1 "qual_file_1"
#define ARG_QUAL_FILE_2 "qual_file_2"
#define ARG_READ_LEN "read_len"
#define ARG_MAX_INDEL "max_indel"
#define ARG_INS_RATE_1 "ins_rate_1"
#define ARG_INS_RATE_2 "ins_rate_2"
#define ARG_DEL_RATE_1 "del_rate_1"
#define ARG_DEL_RATE_2 "del_rate_2"
#define ARG_SEP_FLAG "sep_flag"
#define ARG_PE_FRAG_DIST_MEAN "pe_frag_dist_mean"
#define ARG_PE_FRAG_DIST_STD_DEV "pe_frag_dist_std_dev"
#define ARG_MIN_QUAL "min_qual"
#define ARG_MAX_QUAL "max_qual"
#define ARG_Q_SHIFT_1 "q_shift_1"
#define ARG_Q_SHIFT_2 "q_shift_2"

namespace labw {
namespace art_modern {

    void print_help(const po::options_description& po_desc) { std::cout << po_desc << endl; }

    po::variables_map generate_vm_while_handling_help_version(
        const po::options_description& po_desc, const int argc, char** argv) noexcept
    {
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

        if (vm_.count(ARG_VERSION)) {
            print_version();
            exit(EXIT_SUCCESS);
        }
        if (vm_.count(ARG_HELP)) {
            print_help(po_desc);
            exit(EXIT_SUCCESS);
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
            exit(EXIT_FAILURE);
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
            exit(EXIT_FAILURE);
        }
    }

    INPUT_FILE_TYPE get_input_file_type(const std::string& input_file_type_str, const std::string& input_file_name)
    {
        if (input_file_type_str == INPUT_FILE_TYPE_FASTA) {
            return INPUT_FILE_TYPE::FASTA;
        } else if (input_file_type_str == INPUT_FILE_TYPE_PBSIM3_TEMPLATE) {
            return INPUT_FILE_TYPE::PBSIM3_TEMPLATE;
        } else if (input_file_type_str == INPUT_FILE_TYPE_AUTO) {
            for (const auto& fasta_file_end : std::vector<string> { "fna", "faa", "fa", "fasta" }) {
                if (boost::algorithm::ends_with(input_file_name, fasta_file_end)) {
                    return INPUT_FILE_TYPE::FASTA;
                }
            }
            BOOST_LOG_TRIVIAL(fatal) << "Automatic inference of input file type failed! Modify value of this param (--"
                                     << ARG_INPUT_FILE_TYPE << ") to be one of " << INPUT_FILE_TYPE_FASTA << ", "
                                     << INPUT_FILE_TYPE_PBSIM3_TEMPLATE << ".";
            exit(EXIT_FAILURE);
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "Input file type (--" << ARG_INPUT_FILE_TYPE << ") should be one of "
                                     << INPUT_FILE_TYPE_FASTA << ", " << INPUT_FILE_TYPE_PBSIM3_TEMPLATE << ", "
                                     << INPUT_FILE_TYPE_AUTO << ".";
            exit(EXIT_FAILURE);
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

    INPUT_FILE_PARSER get_input_file_parser(const std::string& input_file_parser_str,
        const std::string& input_file_path, const SIMULATION_MODE simulation_mode)
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
            exit(EXIT_FAILURE);
        }
    }

    std::pair<CoverageInfo, BaseFastaFetch*> get_coverage_info_fasta_fetch(const std::string& fcov_arg_str,
        const INPUT_FILE_TYPE input_file_type, const INPUT_FILE_PARSER input_file_parser,
        const std::string& input_file_name)
    {
        CoverageInfo coverage_info(0);
        BaseFastaFetch* fasta_fetch;
        if (input_file_type == INPUT_FILE_TYPE::PBSIM3_TEMPLATE) {
            // PBSIM3_TEMPLATE have bundled coverage info
        }
        if (fcov_arg_str.empty()) {
            BOOST_LOG_TRIVIAL(fatal) << "Coverage parameter (--" << ARG_FCOV << ") is required.";
            exit(EXIT_FAILURE);
        }

        try {
            auto d = boost::lexical_cast<double>(fcov_arg_str);
            coverage_info = CoverageInfo(d);
        } catch (const boost::bad_lexical_cast&) {
            ifstream X_FOLD(fcov_arg_str, ios::binary);
            coverage_info = CoverageInfo(X_FOLD);
            X_FOLD.close();
        }
        if (input_file_parser == INPUT_FILE_PARSER::STREAM) {
            fasta_fetch = new InMemoryFastaFetch(); // EMpty
        }
        if (input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
            fasta_fetch = new FaidxFetch(input_file_name); // EMpty
        }
        if (input_file_parser == INPUT_FILE_PARSER::MEMORY) {
            if (input_file_type == INPUT_FILE_TYPE::FASTA) {
                fasta_fetch = new InMemoryFastaFetch(input_file_name);
            } else if (input_file_type == INPUT_FILE_TYPE::PBSIM3_TEMPLATE) {
                std::ifstream input_file_stream(input_file_name);
                Pbsim3TranscriptBatcher batcher(std::numeric_limits<int>::max(), input_file_stream);
                fasta_fetch = new InMemoryFastaFetch(batcher.fetch().first);
                coverage_info = batcher.fetch().second;
                input_file_stream.close();
            }
        }
        return { coverage_info, fasta_fetch };
    }

    void shift_emp(map<int, int> map_to_process, const int q_shift, const int min_qual, const int max_qual)
    {
        for (auto& map_to_proces : map_to_process) {
            if (q_shift != 0) {
                if (q_shift < 0 && (-q_shift > map_to_proces.second)) {
                    map_to_proces.second = min_qual;
                } else {
                    map_to_proces.second = min(map_to_proces.second + q_shift, max_qual);
                }
            }
            map_to_proces.second = min(max(map_to_proces.second, min_qual), max_qual);
        }
    }

    void validate_min_max_qual(const int min_qual, const int max_qual)
    {
        if (min_qual < 0 || min_qual > MAX_QUAL) {
            BOOST_LOG_TRIVIAL(fatal) << "Input Error: The minimum quality score must be an integer in [0," << MAX_QUAL
                                     << "]";
            exit(EXIT_FAILURE);
        }
        if (max_qual <= min_qual || max_qual > MAX_QUAL) {
            BOOST_LOG_TRIVIAL(fatal) << "Input Error: The quality score must be an integer in [" << min_qual << ", "
                                     << MAX_QUAL << "]";
            exit(EXIT_FAILURE);
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
            exit(EXIT_FAILURE);
        }
        if (!boost::filesystem::exists(input_file_path)) {
            BOOST_LOG_TRIVIAL(fatal) << "Input file for --" << arg_name << " at '" << input_file_path
                                     << "' does not exist.";
            exit(EXIT_FAILURE);
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
                                         << ", exceeds the maximum first read profile length, " << r1_profile_size
                                         << ".";
            }
            exit(EXIT_FAILURE);
        }

        if ((read_len > r2_profile_size) && art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
            if (r2_profile_size == 0) {
                BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: " << qual_file_2 << ", is not a valid profile.";
            } else {
                BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length, " << read_len
                                         << ", exceeds the maximum second read profile length, " << r2_profile_size
                                         << ".";
            }
            exit(EXIT_FAILURE);
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
        auto seq_file_fai_path = string(fai_path(fasta_path));
        if (!boost::filesystem::exists(boost::filesystem::path(seq_file_fai_path))) {
            BOOST_LOG_TRIVIAL(info) << "Building missing FAI...";
            CExceptionsProxy::requires_numeric(fai_build(fasta_path), USED_HTSLIB_NAME, "Failed to build FAI");
        } else {
            BOOST_LOG_TRIVIAL(info) << "Loading existing FAI...";
            CExceptionsProxy::requires_not_null(
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
            exit(EXIT_FAILURE);
        }
        if (read_len == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Read length must be specified.";
            exit(EXIT_FAILURE);
        }
    }

    void validate_pe_frag_dist(const double pe_frag_dist_mean, const double pe_frag_dist_std_dev, const int read_len,
        const ART_LIB_CONST_MODE art_lib_const_mode, const SIMULATION_MODE art_simulation_mode)
    {
        if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
            if (art_simulation_mode == SIMULATION_MODE::TEMPLATE) {
                if (pe_frag_dist_std_dev != 0 || pe_frag_dist_mean != 0) {
                    BOOST_LOG_TRIVIAL(warning) << "pe_frag_dist_std_dev and "
                                                  "pe_frag_dist_mean ignored for " SIMULATION_MODE_TEMPLATE " mode.";
                }
            } else {
                if (!(pe_frag_dist_std_dev > 0 && pe_frag_dist_mean > 0)) {
                    BOOST_LOG_TRIVIAL(fatal) << "set pe_frag_dist_std_dev and "
                                                "pe_frag_dist_mean for PE reads for " SIMULATION_MODE_WGS
                                                " or " SIMULATION_MODE_TRANS " mode)";
                    exit(EXIT_FAILURE);
                }

                if (pe_frag_dist_mean <= read_len) {
                    BOOST_LOG_TRIVIAL(fatal) << "The read length must be shorter than the "
                                                "pe_frag_dist_mean fragment length specified.";
                    exit(EXIT_FAILURE);
                }
            }
        } else {
            if (art_simulation_mode == SIMULATION_MODE::TEMPLATE
                && (pe_frag_dist_std_dev != 0 || pe_frag_dist_mean != 0)) {
                BOOST_LOG_TRIVIAL(warning) << "pe_frag_dist_std_dev and "
                                              "pe_frag_dist_mean ignored for " ART_LIB_CONST_MODE_SE " mode.";
            }
        }
    }

    void validate_comp_mtx(const INPUT_FILE_PARSER input_file_parser, const SIMULATION_MODE art_simulation_mode,
        const INPUT_FILE_TYPE input_file_type)
    {
        if (input_file_type == INPUT_FILE_TYPE::PBSIM3_TEMPLATE) {
            if (art_simulation_mode == SIMULATION_MODE::WGS) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "Input using PBSim3 transcripts format are not supported for WGS simulation.";
                exit(EXIT_FAILURE);
            } else if (input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
                BOOST_LOG_TRIVIAL(fatal) << "Input using PBSim3 transcripts format are not supported for HTSLib parser";
                exit(EXIT_FAILURE);
            }
        }
        if (art_simulation_mode == SIMULATION_MODE::WGS && input_file_parser == INPUT_FILE_PARSER::STREAM) {
            BOOST_LOG_TRIVIAL(fatal) << "STREAM parser is not supported for WGS simulation ";
            exit(EXIT_FAILURE);
        }
        if (art_simulation_mode != SIMULATION_MODE::WGS && input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
            BOOST_LOG_TRIVIAL(fatal) << "HTSLib parser is only supported in WGS simulation ";
            exit(EXIT_FAILURE);
        }
    }

    ArtParams ArtCmdOpts::parse_args(int argc, char** argv) const
    {
        for (int i = 0; i < argc; i++) {
            args.emplace_back(argv[i]);
        }
        auto vm_ = generate_vm_while_handling_help_version(po_desc_, argc, argv);
        auto art_simulation_mode = get_simulation_mode(vm_[ARG_SIMULATION_MODE].as<string>());
        auto art_lib_const_mode = get_art_lib_const_mode(vm_[ARG_LIB_CONST_MODE].as<string>());
        auto input_file_name = vm_[ARG_INPUT_FILE_NAME].as<string>();
        validate_input_filename(input_file_name, ARG_INPUT_FILE_NAME);
        auto input_file_type = get_input_file_type(vm_[ARG_INPUT_FILE_TYPE].as<string>(), input_file_name);
        auto input_file_parser
            = get_input_file_parser(vm_[ARG_INPUT_FILE_PARSER].as<string>(), input_file_name, art_simulation_mode);
        validate_comp_mtx(input_file_parser, art_simulation_mode, input_file_type);
        if (input_file_parser == INPUT_FILE_PARSER::HTSLIB) {
            validate_htslib_parser(input_file_name);
        }
        auto ci_ff = get_coverage_info_fasta_fetch(
            vm_[ARG_FCOV].as<string>(), input_file_type, input_file_parser, input_file_name);
        auto const& coverage_info = ci_ff.first;
        auto fasta_fetch = ci_ff.second;

        auto id = vm_[ARG_ID].as<string>();

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
        validate_pe_frag_dist(
            pe_frag_dist_mean, pe_frag_dist_std_dev, read_len, art_lib_const_mode, art_simulation_mode);
        auto pe_dist_mean_minus_2_std = static_cast<hts_pos_t>(pe_frag_dist_mean - 2 * pe_frag_dist_std_dev);

        auto parallel = vm_[ARG_PARALLEL].as<int>();
        if (parallel == PARALLEL_ALL) {
            parallel = static_cast<int>(boost::thread::hardware_concurrency());
        } else if (parallel == PARALLEL_DISABLE) {
            parallel = 1;
        }

        auto qdist = read_emp(vm_[ARG_QUAL_FILE_1].as<string>(), vm_[ARG_QUAL_FILE_2].as<string>(), read_len,
            art_lib_const_mode, sep_flag, vm_[ARG_Q_SHIFT_1].as<int>(), vm_[ARG_Q_SHIFT_2].as<int>(),
            vm_[ARG_MIN_QUAL].as<int>(), vm_[ARG_MAX_QUAL].as<int>());
        std::array<double, HIGHEST_QUAL> err_prob {};

        for (int i = 0; i < HIGHEST_QUAL; i++) {
            err_prob[i] = pow(10, -i / 10.0);
        }

        return { art_simulation_mode, art_lib_const_mode, input_file_name, input_file_type, input_file_parser, parallel,
            sep_flag, id, coverage_info, read_len, pe_frag_dist_mean, pe_frag_dist_std_dev, per_base_ins_rate_1,
            per_base_del_rate_1, per_base_ins_rate_2, per_base_del_rate_2, err_prob, pe_dist_mean_minus_2_std, qdist,
            fasta_fetch, out_dispatcher_factory_.create(vm_, fasta_fetch) };
    }

    po::options_description option_parser()
    {
        OutputDispatcherFactory out_dispatcher_factory_;
        out_dispatcher_factory_.add(new FastqReadOutputFactory());
        out_dispatcher_factory_.add(new BamReadOutputFactory());
        out_dispatcher_factory_.add(new HeadlessBamReadOutputFactory());
        po::options_description general_opts("General Options");
        general_opts.add_options()(ARG_HELP, "print out usage information");
        general_opts.add_options()(ARG_VERSION, "display version info");

        po::options_description required_opts("Required Options");
        required_opts.add_options()(ARG_SIMULATION_MODE, po::value<std::string>()->default_value(SIMULATION_MODE_WGS),
            "simulation mode, should be " SIMULATION_MODE_WGS ", " SIMULATION_MODE_TRANS ", " SIMULATION_MODE_TEMPLATE
            ".");
        required_opts.add_options()(ARG_LIB_CONST_MODE, po::value<std::string>()->default_value(ART_LIB_CONST_MODE_SE),
            "library construction mode, should be " ART_LIB_CONST_MODE_SE ", " ART_LIB_CONST_MODE_PE
            ", " ART_LIB_CONST_MODE_MP ".");
        required_opts.add_options()(ARG_INPUT_FILE_PARSER,
            po::value<std::string>()->default_value(INPUT_FILE_PARSER_AUTO),
            "input file parser, should be " INPUT_FILE_PARSER_AUTO ", " INPUT_FILE_PARSER_MEMORY
            ", " INPUT_FILE_PARSER_HTSLIB ", " INPUT_FILE_PARSER_STREAM ".");
        required_opts.add_options()(ARG_INPUT_FILE_TYPE, po::value<std::string>()->default_value(INPUT_FILE_TYPE_AUTO),
            "input file type, should be " INPUT_FILE_TYPE_AUTO ", " INPUT_FILE_TYPE_FASTA
            ", " INPUT_FILE_TYPE_PBSIM3_TEMPLATE ".");

        required_opts.add_options()(ARG_INPUT_FILE_NAME, po::value<std::string>(),
            "the filename of input reference genome, reference "
            "transcriptome, or templates");
        required_opts.add_options()(ARG_FCOV, po::value<std::string>(),
            "the fold of read coverage to be simulated or number of reads/read pairs "
            "generated for each sequence for simulating cDNA reads, or a float for "
            "simulating WGS reads.");

        po::options_description art_opts("ART-specific options");
        art_opts.add_options()(ARG_ID, po::value<std::string>()->default_value(ART_PROGRAM_NAME),
            "the prefix identification tag for read ID");

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
        art_opts.add_options()(ARG_Q_SHIFT_2, po::value<int>()->default_value(0),
            "the amount to shift every second-read quality score by");
        art_opts.add_options()(
            ARG_MIN_QUAL, po::value<int>()->default_value(MIN_QUAL), "the minimum base quality score");
        art_opts.add_options()(
            ARG_MAX_QUAL, po::value<int>()->default_value(MAX_QUAL), "the maximum base quality score");

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

    /*
     *     {
        if (stream) {
            cout << "                  FASTA parser: Stream." << endl;
        } else {
            cout << "                  FASTA parser: HTSLib." << endl;
        }
        if (art_lib_const_mode == ART_LIB_CONST_MODE::SE) {
            cout << "                  Single-end Simulation" << endl;
        } else if (art_lib_const_mode == ART_LIB_CONST_MODE::MP) {
            cout << "                  Matepair-end sequencing simulation" << endl;
        } else {
            cout << "                  Paired-end sequencing simulation" << endl;
        }

        cout << "Parameters used during run" << endl;
        cout << "\tRead Length:\t" << read_len << endl;

        if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
            cout << "\tMean Fragment Length:     " << pe_frag_dist_mean << endl;
            cout << "\tStandard Deviation:       " << pe_frag_dist_std_dev << endl;
        }
        cout << "\tFirst Insertion Rate:     " << ins_rate_1 << endl;
        cout << "\tSecond Insertion Rate:    " << ins_rate_2 << endl;
        cout << "\tFirst Deletion Rate:      " << del_rate_1 << endl;
        cout << "\tSecond Deletion Rate:     " << del_rate_2 << endl;

        cout << "\tFirst quality shift:      " << q_shift_1 << endl;
        cout << "\tSecond quality shift:     " << q_shift_2 << endl;

        if (!sep_flag) {
            cout << "\tProfile Type:             Combined" << endl;
        } else {
            cout << "\tProfile Type:             Separated" << endl;
        }
        cout << "\tID Tag:                   " << id.c_str() << endl << endl;

        cout << "Quality Profile(s)" << endl;

        cout << "\t" << qual_file_1.c_str() << " (user's profile)" << endl;
        if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
            cout << "\t" << qual_file_2.c_str() << " (user's profile)" << endl;
        }
    }
     */
    OutputDispatcherFactory get_output_dispatcher_factory()
    {

        OutputDispatcherFactory out_dispatcher_factory;
        out_dispatcher_factory.add(new FastqReadOutputFactory());
        out_dispatcher_factory.add(new BamReadOutputFactory());
        out_dispatcher_factory.add(new HeadlessBamReadOutputFactory());
        return out_dispatcher_factory;
    }
} // art_modern
} // labw