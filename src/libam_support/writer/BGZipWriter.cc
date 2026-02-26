//
// Created by yuzj on 2/22/26.
//

#include "libam_support/writer/BGZipWriter.hh"

#include "libam_support/utils/mpi_utils.hh"

#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <boost/log/trivial.hpp>

#include <htslib/bgzf.h>

#include <string>

namespace labw::art_modern {
BGZipWriter::BGZipWriter(const std::string& filename, const int level, const std::size_t buffer_size,
    const bool write_gzip, const std::size_t num_threads)
    : WriterInterface(filename)
    , buffer_size_(buffer_size)
{
    gz_file_ = bgzf_open(filename.c_str(), (write_gzip ? "wg" : "w" + std::to_string(level)).c_str());
    if (gz_file_ == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Failed to open file for writing: " << filename;
        abort_mpi();
    }
    bgzf_set_cache_size(gz_file_, static_cast<int>(buffer_size_));
    bgzf_mt(gz_file_, static_cast<int>(num_threads), 0);
}

void BGZipWriter::close()
{
    if (gz_file_ != nullptr) {
        flush();
        const auto errorno = bgzf_close(gz_file_);
        if (errorno != 0) {
            BOOST_LOG_TRIVIAL(error) << "Failed to close BGZF file: " << get_filename() << ", error code: " << errorno;
            abort_mpi();
        }
        gz_file_ = nullptr;
    }
}

std::string BGZipWriter::name() const { return "BGZipWriter"; }

void BGZipWriter::flush()
{
    if (gz_file_ != nullptr) {
        const auto errorno = bgzf_flush(gz_file_);
        if (errorno != 0) {
            BOOST_LOG_TRIVIAL(error) << "Failed to flush BGZF file: " << get_filename() << ", error code: " << errorno;
            abort_mpi();
        }
    }
}
void BGZipWriter::write(std::string&& data)
{
    if (gz_file_ != nullptr) {
        std::string data_copy = std::move(data);
        const auto errorno = bgzf_write(gz_file_, data_copy.data(), data_copy.size());
        if (errorno < 0) {
            BOOST_LOG_TRIVIAL(error) << "Failed to write to BGZF file: " << get_filename()
                                     << ", error code: " << errorno;
            abort_mpi();
        }
    }
}

void BGZipWriter::write(const std::string& data)
{
    if (gz_file_ != nullptr) {
        const auto errorno = bgzf_write(gz_file_, data.data(), data.size());
        if (errorno < 0) {
            BOOST_LOG_TRIVIAL(error) << "Failed to write to BGZF file: " << get_filename()
                                     << ", error code: " << errorno;
            abort_mpi();
        }
    }
}

BGZipWriter::~BGZipWriter() { BGZipWriter::close(); }
} // namespace labw::art_modern
