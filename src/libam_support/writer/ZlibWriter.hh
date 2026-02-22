#if 0
#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include "libam_support/writer/WriterInterface.hh"

#include <zlib.h>

#include <string>

namespace labw::art_modern {
class ZLibWriter final : public WriterInterface {
public:
    explicit ZLibWriter(const std::string& filename, int level, std::size_t buffer_size);

    ~ZLibWriter() override;
    DELETE_COPY(ZLibWriter)
    DELETE_MOVE(ZLibWriter)

    void write(std::string&& data) override;
    void write(const std::string& data) override;
    [[nodiscard]] std::string name() const override;

    void flush() override;

    void close() override;

private:
    gzFile gz_file_;
    std::size_t buffer_size_;
    void update_offset_();
};
} // namespace labw::art_modern
#endif
