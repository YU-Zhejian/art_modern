#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <fstream>

#include <ceu_check/ceu_check_all.h>
#include <ceu_ystrlib/ceu_ystrlib_all.h>
#include <ceu_basic/ceu_c_utils.h>

#include "ArtParams.hh"
#include "art_modern_constants.hh"

using namespace labw::art_modern;
using namespace std;
namespace labw {
namespace art_modern {
    void print_version()
    {
        ceu_ystr_t* full_info = ceu_check_get_full_info();
        char* full_info_cstr = ceu_ystr_to_cstr(full_info);
        cout << full_info_cstr << endl;
        ceu_free_non_null(full_info);
        ceu_free_non_null(full_info_cstr);
    }

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
}
}

void ArtParams::parse_args(const std::vector<std::string>& args)
{
    this->_args = args;
    auto argc = args.size();
    const char* argv[argc];
    int i = 0;
    for (auto const& arg : args) {
        argv[i] = arg.c_str();
        i += 1;
    }

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(static_cast<int>(args.size()), argv, po_desc), vm);
    } catch (...) {
        print_version();
        exit(EXIT_FAILURE);
    }
    po::notify(vm);
    if (vm.count("version")) {
        print_version();
        exit(EXIT_SUCCESS);
    }
    if (vm.count("help")) {
        print_help();
        exit(EXIT_SUCCESS);
    }
    is_amplicon = vm.count("is_amplicon");
    is_mp = vm.count("is_mp");
    no_sam = vm.count("no_sam");
    is_pe = vm.count("is_pe") || is_mp;
    cigar_use_m = vm.count("cigar_use_m");
    sep_flag = vm.count("sep_flag");
    parallel_on_read = vm.count("parallel_on_read");

    seq_file = vm["seq_file"].as<string>();
    out_file_prefix = vm["out_file_prefix"].as<string>();
    read_len = vm["read_len"].as<int>();
    p_cigar.append("read_len").append("=");
    max_indel = vm["max_indel"].as<int>();
    max_indel = vm["max_indel"].as<int>();
    ins_rate_1 = vm["ins_rate_1"].as<double>();
    ins_rate_2 = vm["ins_rate_2"].as<double>();
    del_rate_1 = vm["del_rate_1"].as<double>();
    del_rate_2 = vm["del_rate_2"].as<double>();
    fcov = vm["fcov"].as<string>();
    pe_frag_dist_mean = vm["pe_frag_dist_mean"].as<int>();
    pe_frag_dist_std_dev = vm["pe_frag_dist_std_dev"].as<double>();
    min_qual = vm["min_qual"].as<int>();
    max_qual = vm["max_qual"].as<int>();
    q_shift_1 = vm["q_shift_1"].as<int>();
    q_shift_2 = vm["q_shift_2"].as<int>();
    qual_file_1 = vm["qual_file_1"].as<string>();
    qual_file_2 = vm["qual_file_2"].as<string>();
    id = vm["id"].as<string>();
    max_num_n = vm["max_num_n"].as<int>();
    mask_n = max_num_n > 0;
    parallel = vm["parallel"].as<int>();
}

void ArtParams::validate_args()
{
    if (min_qual < 0 || min_qual > MAX_QUAL) {
        BOOST_LOG_TRIVIAL(fatal) << "Input Error: The minimum quality score must be an integer in [0," << MAX_QUAL << "]" << endl;
        exit(EXIT_FAILURE);
    }
    if (max_qual <= min_qual || max_qual > MAX_QUAL) {
        BOOST_LOG_TRIVIAL(fatal) << "Input Error: The quality score must be an integer in [" << min_qual << ", " << MAX_QUAL << "]" << endl;
        exit(EXIT_FAILURE);
    }
    // Make sure the minimum requirements to run were met if the help tag was not given
    if (fcov.empty()) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: please use only one of the two parameters: fold coverage (-f)." << endl
                                 << endl;
        exit(EXIT_FAILURE);
    }

    if (read_len < 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length must be a positive integer." << endl;
        exit(EXIT_FAILURE);
    }
    if (seq_file.empty() || out_file_prefix.empty() || read_len == 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: An input-file, output-file prefix, read length, and fold coverage must be specified." << endl
                                 << endl;
        exit(EXIT_FAILURE);
    }
    if (boost::filesystem::exists(out_file_prefix)) {
        if (!boost::filesystem::is_directory(out_file_prefix)) {
            BOOST_LOG_TRIVIAL(fatal) << "'" << out_file_prefix << "' is not a directory." << endl
                                     << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        if (!boost::filesystem::create_directories(out_file_prefix)) {
            BOOST_LOG_TRIVIAL(fatal) << "'" << out_file_prefix << "' create directory failed." << endl
                                     << endl;
        }
    }
    if (is_pe && !is_amplicon) {
        if (pe_frag_dist_std_dev > 0 && pe_frag_dist_mean > 0) {
            if (pe_frag_dist_mean >= 2000) {
                is_mp = true;
            } else if (is_mp) {
                BOOST_LOG_TRIVIAL(warning) << "Warning: a mate-pair simulation may be not appropriate for DNA fragment size < 2000bp" << endl;
            }
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "ERROR: set pe_frag_dist_std_dev and pe_frag_dist_mean for PE reads" << endl;
            exit(EXIT_FAILURE);
        }
    } else if (is_amplicon && (pe_frag_dist_std_dev != 0 || pe_frag_dist_mean != 0)) {
        BOOST_LOG_TRIVIAL(warning) << "Warning: pe_frag_dist_std_dev and pe_frag_dist_mean ignored for amplicon." << endl;
    }

    try {
        auto d = boost::lexical_cast<double>(fcov);
        uniform_sequencing_depth = d;
    } catch (const boost::bad_lexical_cast&) {
        ifstream X_FOLD(fcov.c_str(), ios::binary);
        string xfold_line;
        getline(X_FOLD, xfold_line); // Skip 1st line
        while (true) {
            getline(X_FOLD, xfold_line);
            if (xfold_line.empty()) {
                break;
            }
            boost::algorithm::trim(xfold_line);
            auto tab_pos = xfold_line.find('\t');
            sequencing_depth.emplace(xfold_line.substr(0, tab_pos), stod(xfold_line.substr(tab_pos + 1)));
        }
    }

    // Make sure the minimum requirements to run were given for a paired end simulation
    if (is_pe) {
        // use the both default profiles or provide both profiles
        if (qual_file_1.empty() || qual_file_2.empty()) {
            BOOST_LOG_TRIVIAL(fatal) << "Please provide the quality profile of both reads" << endl;
            exit(EXIT_FAILURE);
        }

        if (pe_frag_dist_mean != 0 && pe_frag_dist_std_dev != 0 && pe_frag_dist_mean <= read_len) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length must be shorter than the pe_frag_dist_mean fragment length specified." << endl;
            exit(EXIT_FAILURE);
        }
    }

    if (max_num_n > read_len) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: the cutoff frequency of 'N' for maksing genome region: " << max_num_n << " > the read length ("
                                 << read_len << ")" << endl;
        exit(EXIT_FAILURE);
    }
    per_base_ins_rate_1 = gen_per_base_mutation_rate(read_len, ins_rate_1, max_indel);
    per_base_del_rate_1 = gen_per_base_mutation_rate(read_len, del_rate_1, max_indel);
    per_base_ins_rate_2 = gen_per_base_mutation_rate(read_len, ins_rate_2, max_indel);
    per_base_del_rate_2 = gen_per_base_mutation_rate(read_len, del_rate_2, max_indel);
    for (int i = 0; i < HIGHEST_QUAL; i++) {
        err_prob[i] = pow(10, -i / 10.0);
    }
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

Empdist ArtParams::read_emp() const
{
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
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: " << qual_file_1 << ", is not a valid profile." << endl
                                     << endl;
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length, " << read_len << ", exceeds the maximum first read profile length, "
                                     << r1_profile_size << "." << endl
                                     << endl;
        }
        exit(EXIT_FAILURE);
    }

    if ((read_len > r2_profile_size) && is_pe) {
        if (r2_profile_size == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: " << qual_file_2 << ", is not a valid profile." << endl
                                     << endl;
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: The read length, " << read_len << ", exceeds the maximum second read profile length, "
                                     << r2_profile_size << "." << endl
                                     << endl;
        }
        exit(EXIT_FAILURE);
    }

    if (!sep_flag) {
        for (const auto& i : qdist.qual_dist_first) {
            shift_emp(i, q_shift_1);
        }
        for (const auto& i : qdist.qual_dist_second) {
            shift_emp(i, q_shift_2);
        }
    } else {
        for (const auto& i : qdist.a_qual_dist_first) {
            shift_emp(i, q_shift_1);
        }
        for (const auto& i : qdist.a_qual_dist_second) {
            shift_emp(i, q_shift_2);
        }

        for (const auto& i : qdist.c_qual_dist_first) {
            shift_emp(i, q_shift_1);
        }
        for (const auto& i : qdist.c_qual_dist_second) {
            shift_emp(i, q_shift_2);
        }

        for (const auto& i : qdist.g_qual_dist_first) {
            shift_emp(i, q_shift_1);
        }
        for (const auto& i : qdist.g_qual_dist_second) {
            shift_emp(i, q_shift_2);
        }

        for (const auto& i : qdist.t_qual_dist_first) {
            shift_emp(i, q_shift_1);
        }
        for (const auto& i : qdist.t_qual_dist_second) {
            shift_emp(i, q_shift_2);
        }
    }
    return qdist;
}

string ArtParams::fqfile1(const std::string& gene_name) const
{
    return out_file_prefix + '/' + gene_name + (is_pe ? "_1.fq" : ".fq");
}
string ArtParams::fqfile2(const string& gene_name) const
{
    return is_pe ? out_file_prefix + '/' + gene_name + "_2.fq" : NULL_DEVICE;
}

string ArtParams::samfile(const string& gene_name) const
{
    return out_file_prefix + '/' + gene_name + ".sam";
}

void ArtParams::print_params() const
{
    if (!is_pe) {
        cout << "                  Single-end Simulation" << endl
             << endl;
    } else if (is_mp) {
        cout << "                  Matepair-end sequencing simulation" << endl
             << endl;
    } else {
        cout << "                  Paired-end sequencing simulation" << endl
             << endl;
    }

    cout << "Parameters used during run" << endl;
    cout << "\tRead Length:\t" << read_len << endl;
    if (mask_n) {
        cout << "\tGenome masking 'N' cutoff frequency: \t" << max_num_n << " in " << read_len << endl;
    } else {
        cout << "\t'N' genomic regions masking turned off" << endl;
    }
    if (uniform_sequencing_depth != 0) {
        cout << "\tFold Coverage: (uniform)  " << uniform_sequencing_depth << "X" << endl;
    } else {
        cout << "\tFold Coverage: (file)     " << sequencing_depth.size() << " entries" << endl;
    }

    if (is_pe) {
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
    cout << "\tID Tag:                   " << id.c_str() << endl
         << endl;

    cout << "Quality Profile(s)" << endl;

    cout << "\t" << qual_file_1.c_str() << " (user's profile)" << endl;
    if (is_pe) {
        cout << "\t" << qual_file_2.c_str() << " (user's profile)" << endl;
    }
}

ArtParams::ArtParams()
{
    po_desc.add_options()(
        "qual_file_1", po::value<std::string>()->default_value(""), "the first-read quality profile")(
        "qual_file_2", po::value<std::string>()->default_value(""), "the second-read quality profile. For PE/MP only.")(
        "id", po::value<std::string>()->default_value("ART"), "the prefix identification tag for read ID")(
        "fcov", po::value<std::string>(), "the fold of read coverage to be simulated or number of reads/read pairs generated for each sequence for simulating cDNA reads, or a float for simulating WGS reads.")(
        "help", "print out usage information")(
        "seq_file", po::value<std::string>(), "the filename of input DNA/RNA reference")(
        "ins_rate_1", po::value<double>()->default_value(DEFAULT_INS_RATE_1), "the first-read insertion rate")(
        "ins_rate_2", po::value<double>()->default_value(DEFAULT_INS_RATE_2), "the second-read insertion rate")(
        "del_rate_1", po::value<double>()->default_value(DEFAULT_DEL_RATE_1), "the second-read deletion rate")(
        "del_rate_2", po::value<double>()->default_value(DEFAULT_DEL_RATE_2), "the second-read deletion rate")(
        "sep_flag", "use separate quality profiles for different bases. Default is to use same quality profile regardless its position")(

        "max_indel", po::value<int>()->default_value(DEFAULT_MAX_INDEL), "the maximum total number of insertion and deletion per read")(
        "read_len", po::value<int>(), "read length to be simulated")(
        "pe_frag_dist_mean", po::value<int>()->default_value(0), "Mean distance between DNA/RNA fragments for paired-end simulations")(
        "pe_frag_dist_std_dev", po::value<double>()->default_value(0), "Std. deviation of distance between DNA/RNA fragments for paired-end simulations")(
        "is_mp", "indicate a mate-pair read simulation")(
        "cigar_use_m", "indicate to use CIGAR 'M' instead of '=/X' for alignment match/mismatch")(
        "out_file_prefix", po::value<std::string>()->default_value("art.out"), "the prefix of output filename")(
        "no_sam", "disable syntheis of SAM file")(
        "is_pe", "indicate a paired-end read simulation or to generate reads from both ends. Notice that ART will automatically switch to a mate-pair simulation if the given pe_frag_dist_mean fragment size >= 2000")(
        "q_shift_1", po::value<int>()->default_value(0), "the amount to shift every first-read quality score by")(
        "q_shift_2", po::value<int>()->default_value(0), "the amount to shift every second-read quality score by")(
        "min_qual", po::value<int>()->default_value(MIN_QUAL), "the minimum base quality score")(
        "max_qual", po::value<int>()->default_value(MAX_QUAL), "the maxiumum base quality score")(
        "max_num_n", po::value<int>()->default_value(DEFAULT_MAX_NUM_N), "the cutoff frequency of 'N' in a window size of the read length for masking genomic regions. Use 1 to mask all regions with 'N'. Use 0 to turn off masking")(
        "parallel", po::value<int>()->default_value(PARALLEL_ALL), "Parallel level. -1 for disable, 0 for all CPUs, >=1 to specify number of threads.")("parallel_on_read", "Perform read-level parallelism instead of contig-level")("is_amplicon", "For amplicon simulation")("version", "display version info");
}

void ArtParams::print_help() const
{
    cout << po_desc << endl;
}
