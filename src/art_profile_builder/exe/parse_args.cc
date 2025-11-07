#include "art_profile_builder/exe/parse_args.hh"

#include "art_profile_builder/lib/APBConfig.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Constants.hh"
#include "libam_support/utils/frontend_utils.hh"
#include "libam_support/utils/fs_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/param_utils.hh"
#include "libam_support/utils/seq_utils.hh"

#include <boost/filesystem/operations.hpp>
#include <boost/log/trivial.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <htslib/hts.h>

#include <cstdlib>
#include <string>
#include <vector>

namespace po = boost::program_options;

namespace labw::art_modern {
namespace {
    constexpr char ARG_INPUT_FILE_NAME[] = "i-file";
    constexpr char ARG_PARALLEL[] = "parallel";
    constexpr char ARG_READ_LEN[] = "read_len";
    constexpr char ARG_READ_LEN_1[] = "read_len_1";
    constexpr char ARG_READ_LEN_2[] = "read_len_2";
    constexpr char ARG_NUM_IO_THREADS[] = "i-num_threads";
    constexpr char ARG_IS_PE[] = "is_pe";
    constexpr char ARG_OUT1[] = "o-file1";
    constexpr char ARG_OUT2[] = "o-file2";
    constexpr char ARG_OB[] = "old_behavior";

    htsExactFormat valid_file(const std::string& file_path)
    {
        if (file_path.empty()) {
            BOOST_LOG_TRIVIAL(error) << "Input file path is empty.";
            abort_mpi();
        }
        if (!boost::filesystem::exists(file_path)) {
            BOOST_LOG_TRIVIAL(error) << "Input file does not exist: " << file_path;
            abort_mpi();
        }
        htsFile* file = CExceptionsProxy::assert_not_null(
            hts_open(file_path.c_str(), "r"), USED_HTSLIB_NAME, "Failed to open HTS file.");

        const auto* format = hts_get_format(file);
        if (format->category != sequence_data || format->format == cram || format->format == fasta_format) {
            BOOST_LOG_TRIVIAL(error) << "File is not a sequence data file.";
            hts_close(file);
            abort_mpi();
        }
        hts_close(file);
        return format->format;
    }
    po::options_description option_parser() noexcept
    {
        auto general_opts = general_options();

        po::options_description required_opts("Required Options");
        required_opts.add_options()(
            ARG_INPUT_FILE_NAME, po::value<std::string>(), "the filename of the input FASTQ/SAM/BAM file.");
        required_opts.add_options()(ARG_READ_LEN, po::value<am_read_len_t>(),
            (std::string("maximum read length to be learnt. If the file mode is PE or MP, will use this value on both "
                         "reads. Cannot be specified together with ")
                + ARG_READ_LEN_1 + " or " + ARG_READ_LEN_2)
                .c_str());
        required_opts.add_options()(ARG_READ_LEN_1, po::value<am_read_len_t>(), "read length of read 1 to be learnt");
        required_opts.add_options()(ARG_READ_LEN_2, po::value<am_read_len_t>(), "read length of read 2 to be learnt");

        po::options_description input_flags("Input Flags");
        input_flags.add_options()(ARG_IS_PE, "Whether the input is paired-end. Default: single-end.");
        input_flags.add_options()(ARG_OB,
            "Simulate the behaviour of original ART profile builder. If set, all qualities will be offsetted by 1.");
        input_flags.add_options()(ARG_OUT1, po::value<std::string>(), "Output file name for read 1 profile.");
        input_flags.add_options()(ARG_OUT2, po::value<std::string>(), "Output file name for read 2 profile.");

        po::options_description parallel_opts("Parallelism-related options");
        parallel_opts.add_options()(ARG_PARALLEL, po::value<int>()->default_value(PARALLEL_ALL),
            "Parallel level. -1 for disable, 0 for all CPUs, >=1 to specify number "
            "of threads.");
        const auto io_thread_description
            = std::string("number of threads to use for I/O. Note that every thread specified in ") + ARG_PARALLEL
            + " will create " + ARG_NUM_IO_THREADS + " threads for I/O. ";
        parallel_opts.add_options()(
            ARG_NUM_IO_THREADS, po::value<std::size_t>()->default_value(4), io_thread_description.c_str());

        po::options_description po_desc;
        po_desc.add(general_opts).add(required_opts).add(input_flags).add(parallel_opts);

        return po_desc;
    }
} // namespace

APBConfig parse_args(int argc, char** argv)
{
    auto po_desc_ = option_parser();
    std::vector<std::string> const args { argv, argv + argc };
    BOOST_LOG_TRIVIAL(info) << "ARGS: " << join(args, " ");

    const auto& vm_ = generate_vm_while_handling_help_version(po_desc_, argc, argv);

    const auto is_pe = vm_.count(ARG_IS_PE) > 0;
    const auto is_ob = vm_.count(ARG_OB) > 0;
    am_read_len_t read_len_1 = 0;
    am_read_len_t read_len_2 = 0;
    // Mutal exclusion
    if ((vm_.count(ARG_READ_LEN) > 0) && (vm_.count(ARG_READ_LEN_1) > 0 || vm_.count(ARG_READ_LEN_2) > 0)) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: --" << ARG_READ_LEN << " cannot be specified together with --"
                                 << ARG_READ_LEN_1 << " or --" << ARG_READ_LEN_2 << ".";
        abort_mpi();
    }
    // At least one
    if (vm_.count(ARG_READ_LEN) == 0 && vm_.count(ARG_READ_LEN_1) == 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: --" << ARG_READ_LEN << " or --" << ARG_READ_LEN_1
                                 << " must be specified.";
        abort_mpi();
    }
    // For PE or MP, both read lengths must be specified
    if (is_pe && vm_.count(ARG_READ_LEN_1) > 0 && vm_.count(ARG_READ_LEN_2) == 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: --" << ARG_READ_LEN_2
                                 << " must be specified for PE library construction mode.";
        abort_mpi();
    }
    if (vm_.count(ARG_READ_LEN) != 0) {
        read_len_1 = get_param<am_read_len_t>(vm_, ARG_READ_LEN);
        read_len_2 = is_pe ? read_len_1 : 0;
    } else {
        read_len_1 = get_param<am_read_len_t>(vm_, ARG_READ_LEN_1);
        read_len_2 = is_pe ? get_param<am_read_len_t>(vm_, ARG_READ_LEN_2) : 0;
    }

    const auto input_file_name = get_param<std::string>(vm_, ARG_INPUT_FILE_NAME);
    const auto num_io_threads = get_param<std::size_t>(vm_, ARG_NUM_IO_THREADS);
    const auto num_threads = n_threads_from_parallel(get_param<int>(vm_, ARG_PARALLEL));
    const auto out1 = get_param<std::string>(vm_, ARG_OUT1);
    const auto out2 = is_pe ? get_param<std::string>(vm_, ARG_OUT2) : "";

    const auto format = valid_file(input_file_name);
    prepare_writer(out1);
    if (is_pe) {
        prepare_writer(out2);
    }

    return { input_file_name, read_len_1, read_len_2, num_threads, num_io_threads, is_pe, is_ob, out1, out2, format };
}
} // namespace labw::art_modern
