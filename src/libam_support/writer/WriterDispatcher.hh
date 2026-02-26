#pragma once
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/compression_utils.hh"

#include "libam_support/writer/WriterInterface.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <memory>
#include <string>

namespace labw::art_modern {

class WriterDispatcher {
public:
    DELETE_COPY(WriterDispatcher)
    DELETE_MOVE(WriterDispatcher)
    DEFAULT_DESTRUCTOR(WriterDispatcher)

    static std::unique_ptr<WriterInterface> writer_dispatch(const std::string& filename,
        CompressionType compression_type = CompressionType::GZIP, int compression_level = 6,
        std::size_t buffer_size = 1 << 20, std::size_t num_threads = 1);

    static std::unique_ptr<WriterInterface> parse_args_for_fmt(
        const std::string& fmt_name, const boost::program_options::variables_map& vm, const std::string& resolved_path);
    static void patch_options_for_fmt(const std::string& fmt_name, boost::program_options::options_description& desc);
};
} // namespace labw::art_modern
