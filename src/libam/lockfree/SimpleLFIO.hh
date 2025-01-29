#pragma once
#include "libam/lockfree/LockFreeIO.hh"
#include "libam/utils/class_macros_utils.hh"

#include <memory>
#include <ostream>
#include <string>
#include <utility>

namespace labw::art_modern {

class SimpleLFIO : public LockFreeIO<std::unique_ptr<std::string>> {
public:
    DELETE_MOVE(SimpleLFIO)
    DELETE_COPY(SimpleLFIO)
    void write(std::unique_ptr<std::string> value) override { out_ << *value; }
    explicit SimpleLFIO(std::string name, std::ostream& out)
        : LockFreeIO(std::move(name))
        , out_(out)
    {
    }
    void stop() override
    {
        std::flush(out_);
        LockFreeIO::stop();
    }

private:
    std::ostream& out_;
};
} // namespace labw::art_modern
