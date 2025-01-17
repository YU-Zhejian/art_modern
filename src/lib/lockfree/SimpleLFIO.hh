#pragma once
#include "lockfree/LockFreeIO.hh"

#include <memory>
#include <ostream>
#include <string>

namespace labw::art_modern {

class SimpleLFIO : public LockFreeIO<std::unique_ptr<std::string>> {
public:
    void write(std::unique_ptr<std::string> value) override { out_ << *value; }
    explicit SimpleLFIO(std::ostream& out)
        : out_(out)
    {
    }

private:
    std::ostream& out_;
};
} // namespace labw::art_modern
