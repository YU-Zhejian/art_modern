#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include <iostream>
#include <string>
#include <utility>

namespace labw::art_modern {
class WriterInterface {
public:
    explicit WriterInterface(std::string filename)
        : filename_(std::move(filename))
    {
    }
    virtual ~WriterInterface() = default;
    virtual void write(std::string&& data) = 0;
    virtual void write(const std::string& data) = 0;
    DELETE_COPY(WriterInterface)
    DELETE_MOVE(WriterInterface)
    virtual void close() = 0;
    virtual void flush() = 0;
    [[nodiscard]] std::string get_filename() const { return filename_; }
    [[nodiscard]] virtual std::string name() const = 0;

private:
    std::string filename_;
};
} // namespace labw::art_modern
