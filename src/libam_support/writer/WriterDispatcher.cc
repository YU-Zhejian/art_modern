#include "art_modern_config.h" // NOLINT: for WITH_LIBDEFLATE

#include "libam_support/writer/WriterDispatcher.hh"

#include "libam_support/utils/compression_utils.hh"
#include "libam_support/utils/seq_utils.hh"

#include "libam_support/writer/BGZipWriter.hh"
#include "libam_support/writer/SimpleWriter.hh"
#include "libam_support/writer/WriterInterface.hh"

#if 0
#ifdef WITH_LIBDEFLATE
#include "libam_support/writer/LibDeflateWriter.hh"
#else
#include "libam_support/writer/ZlibWriter.hh"
#endif
#endif

#include <boost/log/trivial.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <fmt/format.h>

#include <memory>
#include <string>

namespace labw::art_modern {
std::unique_ptr<WriterInterface> WriterDispatcher::writer_dispatch(const std::string& filename,
    const CompressionType compression_type, const int compression_level, const std::size_t buffer_size,
    const std::size_t num_threads)
{
    if (compression_type == CompressionType::GZIP) {
#if 0
            if (num_threads > 1)
            {
#endif
        return std::make_unique<BGZipWriter>(filename, compression_level, buffer_size, true, num_threads);
#if 0
            }
#endif
#if 0
#ifdef WITH_LIBDEFLATE
            return std::make_unique<LibDeflateWriter>(filename, compression_level, buffer_size);
#else
            return std::make_unique<ZLibWriter>(filename, compression_level, buffer_size);
#endif
#endif
    }
    if (compression_type == CompressionType::BGZIP) {
        return std::make_unique<BGZipWriter>(filename, compression_level, buffer_size, false, num_threads);
    }
    return std::make_unique<SimpleWriter>(filename, buffer_size);
}

std::unique_ptr<WriterInterface> WriterDispatcher::parse_args_for_fmt(
    const std::string& fmt_name, const boost::program_options::variables_map& vm, const std::string& resolved_path)
{
    const std::string ctype_arg = fmt::format("o-{}-compression", fmt_name);
    const std::string clevel_arg = fmt::format("o-{}-compression_level", fmt_name);
    const std::string bsize_arg = fmt::format("o-{}-buffer_size", fmt_name);
    const std::string nthreads_arg = fmt::format("o-{}-num_threads", fmt_name);
    CompressionType ctype = CompressionType::NOP;
    if (vm.count(ctype_arg) != 0U) {
        const std::string ctype_str = vm[ctype_arg].as<std::string>();
        if (ctype_str == "gzip") {
            ctype = CompressionType::GZIP;
            bool find_gzip_ext = false;
            for (const auto& ext : GZIP_EXTENSIONS) {
                if (ends_with(resolved_path, ext)) {
                    find_gzip_ext = true;
                    break;
                }
            }
            if (!find_gzip_ext) {
                BOOST_LOG_TRIVIAL(warning) << "Specified compression type 'gzip' for format '" << fmt_name
                                           << "' but file extension does not indicate gzip compression. This may lead "
                                              "to issues when downstream tools try to read the file.";
            }
        } else if (ctype_str == "bgzip") {
            ctype = CompressionType::BGZIP;
            bool find_bgzip_ext = false;
            for (const auto& ext : BGZIP_EXTENSIONS) {
                if (ends_with(resolved_path, ext)) {
                    find_bgzip_ext = true;
                    break;
                }
            }
            // Here, GZip extensions should also be considered valud for BGZip.
            for (const auto& ext : GZIP_EXTENSIONS) {
                if (ends_with(resolved_path, ext)) {
                    find_bgzip_ext = true;
                    break;
                }
            }
            if (!find_bgzip_ext) {
                BOOST_LOG_TRIVIAL(warning) << "Specified compression type 'bgzip' for format '" << fmt_name
                                           << "' but file extension does not indicate bgzip compression. This may lead "
                                              "to issues when downstream tools try to read the file.";
            }
        } else if (ctype_str == "none") {
            for (const auto& ext : GZIP_EXTENSIONS) {
                if (ends_with(resolved_path, ext)) {
                    BOOST_LOG_TRIVIAL(warning) << "Specified compression type 'none' for format '" << fmt_name
                                               << "' but file extension indicates gzip compression. This may lead to "
                                                  "issues when downstream tools try to read the file.";
                    break;
                }
            }
            ctype = CompressionType::NONE;
        } else {
            throw std::runtime_error(
                fmt::format("Invalid compression type '{}'. Supported types are 'gzip' and 'none'.", ctype_str));
        }
    } else {
        BOOST_LOG_TRIVIAL(warning) << "Compression type for format '" << fmt_name
                                   << "' not specified. Inferring from file extension.";
        for (const auto& ext : BGZIP_EXTENSIONS) {
            if (ends_with(resolved_path, ext)) {
                ctype = CompressionType::BGZIP;
                break;
            }
        }
        if (ctype == CompressionType::NOP) {
            for (const auto& ext : GZIP_EXTENSIONS) {
                if (ends_with(resolved_path, ext)) {
                    ctype = CompressionType::GZIP;
                    break;
                }
            }
        }
        if (ctype == CompressionType::NOP) {
            ctype = CompressionType::NONE;
        }
    }
    const int clevel = vm.count(clevel_arg) != 0U ? vm[clevel_arg].as<int>() : 6;
    const std::size_t bsize = vm.count(bsize_arg) != 0U ? vm[bsize_arg].as<std::size_t>() : (1 << 20);
    const std::size_t num_threads = vm.count(nthreads_arg) != 0U ? vm[nthreads_arg].as<std::size_t>() : 1;
    return writer_dispatch(resolved_path, ctype, clevel, bsize, num_threads);
}

void WriterDispatcher::patch_options_for_fmt(
    const std::string& fmt_name, boost::program_options::options_description& desc)
{
    const std::string path_arg = fmt::format("o-{}", fmt_name);
    const std::string ctype_arg = fmt::format("o-{}-compression", fmt_name);
    const std::string clevel_arg = fmt::format("o-{}-compression_level", fmt_name);
    const std::string bsize_arg = fmt::format("o-{}-buffer_size", fmt_name);
    const std::string nthreads_arg = fmt::format("o-{}-num_threads", fmt_name);
    desc.add_options()(path_arg.c_str(), boost::program_options::value<std::string>(),
        ("Destination of output " + fmt_name + " file. Unset to disable the writer.").c_str());
    desc.add_options()(ctype_arg.c_str(), boost::program_options::value<std::string>(),
        "Compression type for the output file. Supported values are 'gzip', 'bgzip', and 'none'. If not set, it will "
        "be inferred from the file extension.");
    desc.add_options()(clevel_arg.c_str(), boost::program_options::value<int>()->default_value(6),
        "Compression level for gzip compression. Valid values are typically between 1 (fastest) and 9 (best "
        "compression). Default is 6. Not used when no compression.");
    desc.add_options()(bsize_arg.c_str(), boost::program_options::value<std::size_t>()->default_value(1 << 20),
        "Buffer size in bytes for writing. Default is 1 MiB (1048576 bytes).");
    desc.add_options()(nthreads_arg.c_str(), boost::program_options::value<std::size_t>()->default_value(1),
        "Number of threads to use for compression. Only applicable for gzip compression. Default is 1 (no "
        "multithreading).");
}
} // namespace labw::art_modern
