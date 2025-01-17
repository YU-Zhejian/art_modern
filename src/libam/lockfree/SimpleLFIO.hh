#pragma once
#include "libam/lockfree/LockFreeIO.hh"
#include "libam/utils/class_macros_utils.hh"

#include <memory>
#include <ostream>
#include <string>

namespace labw::art_modern {

class SimpleLFIO : public LockFreeIO<std::unique_ptr<std::string>> {
public:
    DELETE_MOVE(SimpleLFIO)
    DELETE_COPY(SimpleLFIO)
    void write(std::unique_ptr<std::string> value) override { out_ << *value; }
    explicit SimpleLFIO(std::ostream& out)
        : out_(out)
    {
    }
    ~SimpleLFIO() override
    {
        stop();
        std::flush(out_);
    };

private:
    std::ostream& out_;
};
} // namespace labw::art_modern
