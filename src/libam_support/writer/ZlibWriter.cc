#if 0
#include "libam_support/writer/ZlibWriter.hh"

#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#include <zlib.h>

#include <string>

namespace labw::art_modern {
ZLibWriter::~ZLibWriter() { ZLibWriter::close(); }

ZLibWriter::ZLibWriter(const std::string& filename, const int level, const std::size_t buffer_size)
    : WriterInterface(filename)
    , buffer_size_(buffer_size)
{
    gz_file_ = gzopen(filename.c_str(), ("wb" + std::to_string(level)).c_str());
    if (gz_file_ == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Failed to open file for writing: " << filename;
        abort_mpi();
    }
    gzbuffer(gz_file_, static_cast<unsigned int>(buffer_size_));
}

void ZLibWriter::write(std::string&& data)
{
    std::string const data_in = std::move(data);
    const auto bytes_written = gzwrite(gz_file_, data_in.data(), static_cast<unsigned int>(data_in.size()));
    if (bytes_written < 0) {
        BOOST_LOG_TRIVIAL(error) << "Failed to write to gzip file: " << gzerror(gz_file_, nullptr);
        abort_mpi();
    }
    update_offset_();
}

void ZLibWriter::write(const std::string& data)
{
    auto bytes_written = gzwrite(gz_file_, data.data(), static_cast<unsigned int>(data.size()));
    if (bytes_written < 0) {
        BOOST_LOG_TRIVIAL(error) << "Failed to write to gzip file.";
        abort_mpi();
    }
    update_offset_();
}

void ZLibWriter::flush()
{
    auto const errorno = gzflush(gz_file_, Z_SYNC_FLUSH);
    if (errorno != Z_OK) {
        BOOST_LOG_TRIVIAL(error) << "Failed to flush gzip file: " << gzerror(gz_file_, nullptr);
        abort_mpi();
    }
    update_offset_();
}

void ZLibWriter::close()
{
    if (gz_file_ != nullptr) {
        flush();
        const auto errorno = gzclose(gz_file_);
        if (errorno != Z_OK) {
            BOOST_LOG_TRIVIAL(error) << "Failed to close gzip file: " << gzerror(gz_file_, nullptr);
            abort_mpi();
        }
        gz_file_ = nullptr;
    }
}
std::string ZLibWriter::name() const { return "ZLibWriter"; }

void ZLibWriter::update_offset_()
{
    if (gz_file_ == nullptr) {
        return;
    }
    auto current_offset =
#ifdef Z_LARGE64
        gzoffset64
#else
        gzoffset
#endif
        (gz_file_);
    if (current_offset < 0) {
        BOOST_LOG_TRIVIAL(error) << "Failed to write to gzip file.";
        abort_mpi();
    }
}
} // namespace labw::art_modern
#endif
