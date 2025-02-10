#pragma once
#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/fs_utils.hh"

#include <boost/log/trivial.hpp>

#include <atomic>
#include <fstream>
#include <ios>
#include <memory>
#include <ostream>
#include <string>
#include <utility>

namespace labw::art_modern {

class SimpleLFIO : public LockFreeIO<std::unique_ptr<std::string>> {
public:
    DELETE_MOVE(SimpleLFIO)
    DELETE_COPY(SimpleLFIO)
    ~SimpleLFIO() override = default;
    void write(std::unique_ptr<std::string> value) override
    {
        if (closed_) {
            return;
        }
        out_ << *value;
        num_bytes_out_ += value->size();
    }

    SimpleLFIO(std::string name, std::string out_path)
        : LockFreeIO(std::move(name))
        , out_path_(std::move(out_path))
    {
        prepare_writer(out_path_);
        BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: Writer to '" << out_path_ << "' added.";
        out_ = std::ofstream(out_path_, std::ios::out | std::ios::binary);
    }

    SimpleLFIO(std::string name, std::string out_path, const std::string& preamble)
        : SimpleLFIO(std::move(name), std::move(out_path))
    {
        out_ << preamble;
    }



    void flush_and_close() override
    {
        if (closed_) {
            return;
        }
        std::flush(out_);
        out_.close();
        BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: Writer to '" << out_path_ << "' closed.";
        closed_ = true;
    }

private:
    std::atomic<bool> closed_ { false };
    std::ofstream out_;
    std::string out_path_;
};
} // namespace labw::art_modern
