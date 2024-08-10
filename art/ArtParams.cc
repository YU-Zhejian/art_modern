#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <fstream>

#include "ArtParams.hh"
#include "art_modern_constants.hh"

using namespace labw::art_modern;
using namespace std;
namespace labw {
namespace art_modern {
    void print_version()
    {
        // TODO
    }

    std::vector<double> gen_per_base_mutation_rate(int read_len, double p,
        int max_num)
    {
        std::vector<double> rate;
        if (max_num == 0 || p < 1E-30) {
            return rate;
        }

        double tp;
        double p_cdf = 0;
        for (size_t i = 0; i < read_len; i++) {
            tp = boost::math::cdf(
                boost::math::complement(boost::math::binomial(read_len, p), i));
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
} // namespace art_modern
} // namespace labw

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

    po::variables_map vm;
    try {
        po::store(
            po::parse_command_line(static_cast<int>(args.size()), argv, po_desc),
            vm);
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

    auto simulation_mode_str = vm["mode"].as<string>();
    if (simulation_mode_str == "wgs") {
        art_simulation_mode = ART_SIMULATION_MODE::WGS;
    } else if (simulation_mode_str == "trans") {
        art_simulation_mode = ART_SIMULATION_MODE::TRANS;
    } else if (simulation_mode_str == "template") {
        art_simulation_mode = ART_SIMULATION_MODE::TEMPLATE;
    } else {
        BOOST_LOG_TRIVIAL(fatal)
            << R"(Simulation mode (--mode) should be one of "wgs", "trans" and "template")";
        exit(EXIT_FAILURE);
    }
    auto lib_const_mode_str = vm["lc"].as<string>();
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
    no_sam = vm.count("no_sam");
    stream = vm.count("stream");
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
        BOOST_LOG_TRIVIAL(fatal)
            << "Input Error: The minimum quality score must be an integer in [0,"
            << MAX_QUAL << "]";
        exit(EXIT_FAILURE);
    }
    if (max_qual <= min_qual || max_qual > MAX_QUAL) {
        BOOST_LOG_TRIVIAL(fatal)
            << "Input Error: The quality score must be an integer in [" << min_qual
            << ", " << MAX_QUAL << "]";
        exit(EXIT_FAILURE);
    }
    // Make sure the minimum requirements to run were met if the help tag was not
    // given
    if (fcov.empty()) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: please use only one of the two "
                                    "parameters: fold coverage (-f).";
        exit(EXIT_FAILURE);
    }

    if (read_len < 0) {
        BOOST_LOG_TRIVIAL(fatal)
            << "Fatal Error: The read length must be a positive integer.";
        exit(EXIT_FAILURE);
    }
    if (seq_file.empty() || out_file_prefix.empty() || read_len == 0) {
        BOOST_LOG_TRIVIAL(fatal)
            << "Fatal Error: An input-file, output-file prefix, read length, and "
               "fold coverage must be specified.";
        exit(EXIT_FAILURE);
    }
    if (boost::filesystem::exists(out_file_prefix)) {
        if (!boost::filesystem::is_directory(out_file_prefix)) {
            BOOST_LOG_TRIVIAL(fatal)
                << "'" << out_file_prefix << "' is not a directory.";
            exit(EXIT_FAILURE);
        }
    } else {
        if (!boost::filesystem::create_directories(out_file_prefix)) {
            BOOST_LOG_TRIVIAL(fatal)
                << "'" << out_file_prefix << "' create directory failed.";
            exit(EXIT_FAILURE);
        }
    }
    if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        if (art_simulation_mode != ART_SIMULATION_MODE::TEMPLATE && !(pe_frag_dist_std_dev > 0 && pe_frag_dist_mean > 0)) {
            BOOST_LOG_TRIVIAL(fatal) << "ERROR: set pe_frag_dist_std_dev and "
                                        "pe_frag_dist_mean for PE reads for \"wgs\" or \"trans\" mode";
            exit(EXIT_FAILURE);
        }
        if (art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE && (pe_frag_dist_std_dev != 0 || pe_frag_dist_mean != 0)) {
            BOOST_LOG_TRIVIAL(warning) << "Warning: pe_frag_dist_std_dev and "
                                          "pe_frag_dist_mean ignored for template mode.";
        }
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
            sequencing_depth.emplace(xfold_line.substr(0, tab_pos),
                stod(xfold_line.substr(tab_pos + 1)));
        }
    }

    // Make sure the minimum requirements to run were given for a paired end
    // simulation
    if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        // use the both default profiles or provide both profiles
        if (qual_file_1.empty() || qual_file_2.empty()) {
            BOOST_LOG_TRIVIAL(fatal)
                << "Please provide the quality profile of both reads";
            exit(EXIT_FAILURE);
        }

        if (pe_frag_dist_mean != 0 && pe_frag_dist_std_dev != 0 && pe_frag_dist_mean <= read_len) {
            BOOST_LOG_TRIVIAL(fatal)
                << "Fatal Error: The read length must be shorter than the "
                   "pe_frag_dist_mean fragment length specified.";
            exit(EXIT_FAILURE);
        }
    }

    if (max_num_n > read_len) {
        BOOST_LOG_TRIVIAL(fatal)
            << "Error: the cutoff frequency of 'N' for maksing genome region: "
            << max_num_n << " > the read length (" << read_len << ")";
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
            BOOST_LOG_TRIVIAL(fatal)
                << "Fatal Error: " << qual_file_1 << ", is not a valid profile.";
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << "Fatal Error: The read length, " << read_len
                << ", exceeds the maximum first read profile length, "
                << r1_profile_size << ".";
        }
        exit(EXIT_FAILURE);
    }

    if ((read_len > r2_profile_size) && art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        if (r2_profile_size == 0) {
            BOOST_LOG_TRIVIAL(fatal)
                << "Fatal Error: " << qual_file_2 << ", is not a valid profile.";
        } else {
            BOOST_LOG_TRIVIAL(fatal)
                << "Fatal Error: The read length, " << read_len
                << ", exceeds the maximum second read profile length, "
                << r2_profile_size << ".";
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
    return out_file_prefix + '/' + gene_name + (art_lib_const_mode != ART_LIB_CONST_MODE::SE ? "_1.fq" : ".fq");
}
string ArtParams::fqfile2(const string& gene_name) const
{
    return art_lib_const_mode != ART_LIB_CONST_MODE::SE ? out_file_prefix + '/' + gene_name + "_2.fq" : "";
}

string ArtParams::samfile(const string& gene_name) const
{
    return out_file_prefix + '/' + gene_name + ".sam";
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
    if (mask_n) {
        cout << "\tGenome masking 'N' cutoff frequency: \t" << max_num_n << " in "
             << read_len << endl;
    } else {
        cout << "\t'N' genomic regions masking turned off" << endl;
    }
    if (uniform_sequencing_depth != 0) {
        cout << "\tFold Coverage: (uniform)  " << uniform_sequencing_depth << "X"
             << endl;
    } else {
        cout << "\tFold Coverage: (file)     " << sequencing_depth.size()
             << " entries" << endl;
    }

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
    cout << "\tID Tag:                   " << id.c_str() << endl
         << endl;

    cout << "Quality Profile(s)" << endl;

    cout << "\t" << qual_file_1.c_str() << " (user's profile)" << endl;
    if (art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
        cout << "\t" << qual_file_2.c_str() << " (user's profile)" << endl;
    }
}

ArtParams::ArtParams()
{
    po_desc.add_options()("help", "print out usage information");
    po_desc.add_options()("version", "display version info");

    po_desc.add_options()(
        "mode", po::value<std::string>()->default_value("wgs"),
        R"(simulation mode, should be "wgs", "trans" or "template")");
    po_desc.add_options()(
        "lc", po::value<std::string>()->default_value("se"),
        R"(library construction mode, should be "se", "pe" or "mp")");

    po_desc.add_options()("seq_file", po::value<std::string>(),
        "the filename of input reference genome, reference "
        "transcriptome, or templates");
    po_desc.add_options()(
        "fcov", po::value<std::string>(),
        "the fold of read coverage to be simulated or number of reads/read pairs "
        "generated for each sequence for simulating cDNA reads, or a float for "
        "simulating WGS reads.");
    po_desc.add_options()("out_file_prefix",
        po::value<std::string>()->default_value("art.out"),
        "the prefix of output filename");

    po_desc.add_options()("id", po::value<std::string>()->default_value("ART"),
        "the prefix identification tag for read ID");

    po_desc.add_options()("qual_file_1",
        po::value<std::string>()->default_value(""),
        "the first-read quality profile");
    po_desc.add_options()("qual_file_2",
        po::value<std::string>()->default_value(""),
        "the second-read quality profile. For PE/MP only.");
    po_desc.add_options()("ins_rate_1",
        po::value<double>()->default_value(DEFAULT_INS_RATE_1),
        "the first-read insertion rate");
    po_desc.add_options()("ins_rate_2",
        po::value<double>()->default_value(DEFAULT_INS_RATE_2),
        "the second-read insertion rate");
    po_desc.add_options()("del_rate_1",
        po::value<double>()->default_value(DEFAULT_DEL_RATE_1),
        "the second-read deletion rate");
    po_desc.add_options()("del_rate_2",
        po::value<double>()->default_value(DEFAULT_DEL_RATE_2),
        "the second-read deletion rate");
    po_desc.add_options()(
        "sep_flag", "use separate quality profiles for different bases. Default "
                    "is to use same quality profile regardless its position");
    po_desc.add_options()(
        "max_indel", po::value<int>()->default_value(DEFAULT_MAX_INDEL),
        "the maximum total number of insertion and deletion per read");
    po_desc.add_options()("read_len", po::value<int>(),
        "read length to be simulated");
    po_desc.add_options()(
        "pe_frag_dist_mean", po::value<int>()->default_value(0),
        "Mean distance between DNA/RNA fragments for paired-end simulations");
    po_desc.add_options()(
        "pe_frag_dist_std_dev", po::value<double>()->default_value(0),
        "Std. deviation of distance between DNA/RNA fragments for paired-end "
        "simulations");
    po_desc.add_options()(
        "cigar_use_m", "indicate to use CIGAR 'M' instead of '=/X' for alignment "
                       "match/mismatch");
    po_desc.add_options()("no_sam", "disable syntheis of SAM file");
    po_desc.add_options()(
        "q_shift_1", po::value<int>()->default_value(0),
        "the amount to shift every first-read quality score by");
    po_desc.add_options()(
        "q_shift_2", po::value<int>()->default_value(0),
        "the amount to shift every second-read quality score by");
    po_desc.add_options()("min_qual", po::value<int>()->default_value(MIN_QUAL),
        "the minimum base quality score");
    po_desc.add_options()("max_qual", po::value<int>()->default_value(MAX_QUAL),
        "the maxiumum base quality score");
    po_desc.add_options()(
        "max_num_n", po::value<int>()->default_value(DEFAULT_MAX_NUM_N),
        "the cutoff frequency of 'N' in a window size of the read length for "
        "masking genomic regions. Use 1 to mask all regions with 'N'. Use 0 to "
        "turn off masking");
    po_desc.add_options()(
        "parallel", po::value<int>()->default_value(PARALLEL_ALL),
        "Parallel level. -1 for disable, 0 for all CPUs, >=1 to specify number "
        "of threads.");
    po_desc.add_options()(
        "parallel_on_read",
        "Perform read-level parallelism instead of contig-level");
    po_desc.add_options()(
        "stream", "If specified, will use streamline FASTA parser (For "
                  "transcriptome/amplicon); Otherwise will use HTSLib indexed "
                  "FASTA parser (For genome).");
}

void ArtParams::print_help() const { cout << po_desc << endl; }
