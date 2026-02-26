#if 0
#pragma once

#include "art_modern_config.h" // NOLINT: for WITH_LIBDEFLATE

#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include "libam_support/writer/WriterInterface.hh"

#include <libdeflate.h>

#include <fstream>
#include <string>
#include <vector>

namespace labw::art_modern {
#ifdef WITH_LIBDEFLATE
class LibDeflateWriter final : public WriterInterface {
public:
    // Compression levels usually range from 0 to 12
    explicit LibDeflateWriter(const std::string& filename, int level, std::size_t buffer_size);

    ~LibDeflateWriter() override;

    DELETE_COPY(LibDeflateWriter)
    DELETE_MOVE(LibDeflateWriter)

    // Write moves the string into our internal buffer
    void write(std::string&& data) override;
    void write(const std::string& data) override;

    // Manual flush: Compresses whatever is currently in the buffer
    void flush() override;

    void close() override;

    [[nodiscard]] std::string name() const override;

private:
    void flush_buffer_to_disk();

    std::ofstream out_stream_;
    std::vector<char> buffer_;
    libdeflate_compressor* compressor_;
    std::size_t compressed_size_ = 0;
    std::vector<char> compressed_out_;
    std::size_t buffer_size_;
};
#endif

} // namespace labw::art_modern
#endif
