#if 0
#include "libam_support/writer/LibDeflateWriter.hh"

#include "art_modern_config.h" // NOLINT: for WITH_LIBDEFLATE

#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#include <libdeflate.h>

#include <fstream>
#include <string>
#include <vector>

namespace labw::art_modern
{
#ifdef WITH_LIBDEFLATE
    LibDeflateWriter::LibDeflateWriter(const std::string& filename, const int level, const std::size_t buffer_size)
        : WriterInterface(filename),
    compressor_(libdeflate_alloc_compressor(level))
        , buffer_size_(buffer_size)
    {
        out_stream_.open(filename, std::ios::out | std::ios::binary);
        if (!out_stream_) {
            BOOST_LOG_TRIVIAL(error) << "Failed to open file for writing: " << filename;
            abort_mpi();
        }
        buffer_.reserve(buffer_size_);
        if (compressor_ == nullptr) {
            BOOST_LOG_TRIVIAL(error) << "Failed to allocate compressor with level " << level;
            abort_mpi();
        }
    }

    LibDeflateWriter::~LibDeflateWriter()
    {
        LibDeflateWriter::close();
    }

    void LibDeflateWriter::write(std::string&& data)
    {
        const std::string data_in = std::move(data);
        buffer_.insert(buffer_.end(), data_in.begin(), data_in.end());
        if (buffer_.size() >= buffer_size_) {
            flush_buffer_to_disk();
        }
    }

    void LibDeflateWriter::write(const std::string& data)
    {
        buffer_.insert(buffer_.end(), data.begin(), data.end());
        if (buffer_.size() >= buffer_size_) {
            flush_buffer_to_disk();
        }
    }

    void LibDeflateWriter::flush()
    {
        if (!buffer_.empty()) {
            flush_buffer_to_disk();
        }
        out_stream_.flush();
    }

    void LibDeflateWriter::close()
    {
        if (out_stream_.is_open()) {
            flush();
            out_stream_.close();
            libdeflate_free_compressor(compressor_);
        }
    } std::string LibDeflateWriter::name() const { return "LibDeflateWriter"; }


    void LibDeflateWriter::flush_buffer_to_disk()
    {
        if (buffer_.empty()) {
            return;
        }

        // libdeflate needs a destination buffer.
        // Worst case: compressed size slightly larger than input
        const size_t max_compressed_size = libdeflate_gzip_compress_bound(compressor_, buffer_.size());
        compressed_out_.resize(max_compressed_size);

        const size_t actual_size = libdeflate_gzip_compress(
            compressor_, buffer_.data(), buffer_.size(), compressed_out_.data(), compressed_out_.size());

        if (actual_size > 0) {
            // Write the compressed size first (common for block-based custom formats)
            // or just write the raw block if you're using a specific container format.
            out_stream_.write(compressed_out_.data(), static_cast<std::streamsize>(actual_size));
        } else {
            BOOST_LOG_TRIVIAL(error) << "Compression failed for the current buffer.";
            abort_mpi();
        }
        compressed_size_ += actual_size;
        buffer_.clear();
    }
#endif
} // namespace labw::art_modern
#endif
