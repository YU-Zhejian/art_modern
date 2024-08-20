#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <fstream>
#include "ArtCmdOpts.hh"

#include <htslib/hfile.h>
#include <htslib/hts.h>

#include "ArtParams.hh"
#include "art_modern_constants.hh"
#include "global_variables.hh"
#include "out/BamReadOutput.hh"
#include "out/HeadlessBamReadOutput.hh"
#include "out/FastqReadOutput.hh"
#include "CExceptionsProxy.hh"

#include <fasta/InMemoryFastaFetch.hh>

using namespace std;
namespace po = boost::program_options;

namespace labw {
namespace art_modern {


    std::vector<double> gen_per_base_mutation_rate(int read_len, double p, int max_num)
    {
        std::vector<double> rate;
        if (max_num == 0 || p < 1E-30) {
            return rate;
        }

        double tp;
        double p_cdf = 0;
        for (size_t i = 0; i < read_len; i++) {
            tp = boost::math::cdf(boost::math::complement(boost::math::binomial(read_len, p), i));
            rate.push_back(tp);
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

    void ArtParams::parse_args(const std::vector<std::string>& args)
    {
        this->_args = args;
        auto argc = args.size();
        BOOST_LOG_TRIVIAL(info) << "ARGS: " << boost::algorithm::join(args, " ");
        if (argc <= 1) {
            print_help();
            exit(EXIT_FAILURE);
        }
        const char* argv[argc];
        int i = 0;
        for (auto const& arg : args) {
            argv[i] = arg.c_str();
            i += 1;
        }

        try {
            po::store(po::parse_command_line(static_cast<int>(args.size()), argv, po_desc), vm_);
        } catch (const exception& exp) {
            BOOST_LOG_TRIVIAL(fatal) << exp.what();
            print_help();
            exit(EXIT_FAILURE);
        }
        po::notify(vm_);

        if (vm_.count("version")) {
            print_version();
            exit(EXIT_SUCCESS);
        }
        if (vm_.count("help")) {
            print_help();
            exit(EXIT_SUCCESS);
        }

        auto simulation_mode_str = vm_["mode"].as<string>();
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
        if (lib_const_mode_str == "se") {
            art_lib_const_mode = ART_LIB_CONST_MODE::SE;
        } else if (lib_const_mode_str == "pe") {
            art_lib_const_mode = ART_LIB_CONST_MODE::PE;
        } else if (lib_const_mode_str == "mp") {
            art_lib_const_mode = ART_LIB_CONST_MODE::MP;
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << R"(Library construction mode (--lc) should be one of "wgs", "trans" and "template")";
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

    void validate_min_max_qual(const ArtParams& art_params)  {
        if (art_params.min_qual < 0 || art_params.min_qual > MAX_QUAL) {
            BOOST_LOG_TRIVIAL(fatal) << "Input Error: The minimum quality score must be an integer in [0," << MAX_QUAL
                                     << "]";
            exit(EXIT_FAILURE);
        }
        if (art_params.max_qual <= art_params.min_qual || art_params.max_qual > MAX_QUAL) {
            BOOST_LOG_TRIVIAL(fatal) << "Input Error: The quality score must be an integer in [" << art_params.min_qual << ", "
                                     << MAX_QUAL << "]";
            exit(EXIT_FAILURE);
        }
    }

void validate_fasta_parser(const ArtParams& art_params){
    if (!art_params.stream) {
        const char* fasta_path = art_params.seq_file.c_str();
        BOOST_LOG_TRIVIAL(info) << "HTSLib parser requested. Checking FAI...";
        auto seq_file_fai_path = string(fai_path(fasta_path));
        if (!boost::filesystem::exists(boost::filesystem::path(seq_file_fai_path))) {
            BOOST_LOG_TRIVIAL(info) << "Building missing FAI...";
            CExceptionsProxy::requires_numeric(
                fai_build(fasta_path), USED_HTSLIB_NAME, "Failed to build FAI");
        } else {
            BOOST_LOG_TRIVIAL(info) << "Loading existing FAI...";
            CExceptionsProxy::requires_not_null(
                fai_load_format(fasta_path, FAI_FASTA), USED_HTSLIB_NAME, "Failed to load FAI");
        }
    }
    }

    void ArtParams::validate_args()
    {
        validate_min_max_qual(*this);
        if (read_len < 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length must be a positive integer.";
            exit(EXIT_FAILURE);
        }
        if (seq_file.empty() || read_len == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: An input-file, read length must be specified.";
            exit(EXIT_FAILURE);
        }
        if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
            if (art_simulation_mode != SIMULATION_MODE::TEMPLATE
                && !(pe_frag_dist_std_dev > 0 && pe_frag_dist_mean > 0)) {
                BOOST_LOG_TRIVIAL(fatal) << R"(set pe_frag_dist_std_dev and "
                                            "pe_frag_dist_mean for PE reads for "wgs" or "trans" mode)";
                exit(EXIT_FAILURE);
            }
            if (art_simulation_mode == SIMULATION_MODE::TEMPLATE
                && (pe_frag_dist_std_dev != 0 || pe_frag_dist_mean != 0)) {
                BOOST_LOG_TRIVIAL(warning) << R"(pe_frag_dist_std_dev and "
                                              "pe_frag_dist_mean ignored for "template" mode.)";
            }
        }


        // Make sure the minimum requirements to run were given for a paired end
        // simulation
        if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
            // use the both default profiles or provide both profiles
            if (qual_file_1.empty() || qual_file_2.empty()) {
                BOOST_LOG_TRIVIAL(fatal) << "Please provide the quality profile of both reads";
                exit(EXIT_FAILURE);
            }

            if (pe_frag_dist_mean != 0 && pe_frag_dist_std_dev != 0 && pe_frag_dist_mean <= read_len) {
                BOOST_LOG_TRIVIAL(fatal) << "The read length must be shorter than the "
                                            "pe_frag_dist_mean fragment length specified.";
                exit(EXIT_FAILURE);
            }
        }

        per_base_ins_rate_1 = gen_per_base_mutation_rate(read_len, ins_rate_1, max_indel);
        per_base_del_rate_1 = gen_per_base_mutation_rate(read_len, del_rate_1, max_indel);
        per_base_ins_rate_2 = gen_per_base_mutation_rate(read_len, ins_rate_2, max_indel);
        per_base_del_rate_2 = gen_per_base_mutation_rate(read_len, del_rate_2, max_indel);

        for (int i = 0; i < HIGHEST_QUAL; i++) {
            err_prob[i] = pow(10, -i / 10.0);
        }

        pe_dist_mean_minus_2_std = static_cast<hts_pos_t>(pe_frag_dist_mean - 2 * pe_frag_dist_std_dev);

        validate_fasta_parser(*this);
        read_emp();
    }

    void ArtParams::shift_emp(map<int, int> map_to_process, int q_shift) const
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

    void ArtParams::read_emp()
    {
        qdist_ = Empdist(qual_file_1, qual_file_2, sep_flag);
        size_t r1_profile_size;
        size_t r2_profile_size;
        if (sep_flag) {
            r1_profile_size = qdist_.a_qual_dist_first.size();
            r2_profile_size = qdist_.a_qual_dist_second.size();
        } else {
            r1_profile_size = qdist_.qual_dist_first.size();
            r2_profile_size = qdist_.qual_dist_second.size();
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

        if (!sep_flag) {
            for (const auto& i : qdist_.qual_dist_first) {
                shift_emp(i, q_shift_1);
            }
            for (const auto& i : qdist_.qual_dist_second) {
                shift_emp(i, q_shift_2);
            }
        } else {
            for (const auto& i : qdist_.a_qual_dist_first) {
                shift_emp(i, q_shift_1);
            }
            for (const auto& i : qdist_.a_qual_dist_second) {
                shift_emp(i, q_shift_2);
            }

            for (const auto& i : qdist_.c_qual_dist_first) {
                shift_emp(i, q_shift_1);
            }
            for (const auto& i : qdist_.c_qual_dist_second) {
                shift_emp(i, q_shift_2);
            }

            for (const auto& i : qdist_.g_qual_dist_first) {
                shift_emp(i, q_shift_1);
            }
            for (const auto& i : qdist_.g_qual_dist_second) {
                shift_emp(i, q_shift_2);
            }

            for (const auto& i : qdist_.t_qual_dist_first) {
                shift_emp(i, q_shift_1);
            }
            for (const auto& i : qdist_.t_qual_dist_second) {
                shift_emp(i, q_shift_2);
            }
        }
    }

    void ArtParams::print_params() const
    {
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

    ArtParams::ArtParams(): coverage_info(0.0)
    {
        return;
    }

    void ArtParams::print_help() const { po_desc.print(cout, 0); }
    std::shared_ptr<BaseReadOutput> ArtParams::get_output_dispatcher() const
    {
        // FIXME:
        auto ff = std::make_shared<InMemoryFastaFetch>();
        return out_dispatcher_factory_.create(vm_, ff);
    }

void ArtParams::read_fcov()
{
    if (fcov.empty()) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: please use only one of the two "
                                    "parameters: fold coverage (-f).";
        exit(EXIT_FAILURE);
    }

    try {
        auto d = boost::lexical_cast<double>(fcov);
        coverage_info = CoverageInfo(d);
    } catch (const boost::bad_lexical_cast&) {
        ifstream X_FOLD(fcov.c_str(), ios::binary);
        coverage_info = CoverageInfo(X_FOLD);
        X_FOLD.close();
    }
}

} // namespace art_modern
} // namespace labw
