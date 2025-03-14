#include "art_modern_config.h" // NOLINT: For WITH_MPI

#include "libam_support/utils/fs_utils.hh"

#include "libam_support/Dtypes.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/log/trivial.hpp>

#include <string>

namespace labw::art_modern {
void validate_input_filename(const std::string& input_file_path, const std::string& arg_name)
{
    if (input_file_path.empty()) {
        BOOST_LOG_TRIVIAL(fatal) << "An input file path for --" << arg_name << " must be specified.";
        abort_mpi();
    }
    if (!boost::filesystem::exists(input_file_path)) {
        BOOST_LOG_TRIVIAL(fatal) << "Input file for --" << arg_name << " at '" << input_file_path
                                 << "' does not exist.";
        abort_mpi();
    }
    if (!boost::filesystem::is_regular_file(input_file_path)) {
        BOOST_LOG_TRIVIAL(warning) << "Input file for --" << arg_name << " at '" << input_file_path
                                   << "' is not a regular file.";
#ifdef WITH_MPI
        BOOST_LOG_TRIVIAL(fatal) << "Irregular file is NOT allowed under MPI.";
        abort_mpi();
#endif
    }
}
am_filelen_t get_file_size(const std::string& file_path) noexcept
{
    if (!boost::filesystem::is_regular_file(file_path)) {
        return -1;
    }

    try {
        return static_cast<am_filelen_t>(boost::filesystem::file_size(file_path));
    } catch (const boost::filesystem::filesystem_error&) {
        return -1;
    }
}
void ensure_directory_exists(const std::string& dir_path)
{
    if (boost::filesystem::exists(dir_path)) {
        if (boost::filesystem::is_directory(dir_path)) {
            return;
        }

        BOOST_LOG_TRIVIAL(fatal) << "Path '" << dir_path << "' exists but is not a directory.";
        abort_mpi();
    }
    BOOST_LOG_TRIVIAL(info) << "Boost::filesystem: mkdir -p '" << dir_path << "'";
    boost::filesystem::create_directories(dir_path);
}
void prepare_writer(const std::string& output_file_path)
{
    ensure_directory_exists(boost::filesystem::path(output_file_path).parent_path().string());
}
} // namespace labw::art_modern
// labw