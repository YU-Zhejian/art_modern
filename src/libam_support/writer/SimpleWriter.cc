#include "libam_support/writer/SimpleWriter.hh"

#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#include <fstream>
#include <string>

namespace labw::art_modern {
SimpleWriter::SimpleWriter(const std::string& filename, const std::size_t buffer_size)
    : WriterInterface(filename)
    , buffer_(buffer_size)
{

    ofs_.rdbuf()->pubsetbuf(buffer_.data(), static_cast<long>(buffer_size));
    ofs_.open(filename, std::ios::out | std::ios::binary);

    if (!ofs_) {
        BOOST_LOG_TRIVIAL(error) << "Failed to open file for writing: " << filename;
        abort_mpi();
    }
}
std::string SimpleWriter::name() const { return "SimpleWriter"; }

SimpleWriter::~SimpleWriter() { SimpleWriter::close(); }

void SimpleWriter::write(std::string&& data)
{
    std::string const data_in = std::move(data);
    ofs_ << data_in;
}

void SimpleWriter::close()
{
    if (ofs_.is_open()) {
        ofs_.close();
    }
}

void SimpleWriter::flush() { ofs_.flush(); }

void SimpleWriter::write(const std::string& data) { ofs_ << data; }
} // namespace labw::art_modern
