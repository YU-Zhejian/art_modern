#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include "libam_support/writer/WriterInterface.hh"

#include <htslib/bgzf.h>

#include <string>

namespace labw::art_modern {
class BGZipWriter final : public WriterInterface {
public:
    explicit BGZipWriter(
        const std::string& filename, int level, std::size_t buffer_size, bool write_gzip, std::size_t num_threads);

    ~BGZipWriter() override;
    DELETE_COPY(BGZipWriter)
    DELETE_MOVE(BGZipWriter)

    void write(std::string&& data) override;
    void write(const std::string& data) override;
    void flush() override;
    void close() override;
    [[nodiscard]] std::string name() const override;

private:
    BGZF* gz_file_;
    std::size_t buffer_size_;
};
} // namespace labw::art_modern
