#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include "libam_support/writer/WriterInterface.hh"

#include <fstream>
#include <string>
#include <vector>

namespace labw::art_modern {
class SimpleWriter : public WriterInterface {
public:
    DELETE_COPY(SimpleWriter)
    DELETE_MOVE(SimpleWriter)
    explicit SimpleWriter(const std::string& filename, std::size_t buffer_size);
    ~SimpleWriter() override;

    void write(std::string&& data) override;
    void close() override;
    void flush() override;
    void write(const std::string& data) override;
    [[nodiscard]] std::string name() const override;

private:
    std::ofstream ofs_;
    std::vector<char> buffer_;
};
} // namespace labw::art_modern
