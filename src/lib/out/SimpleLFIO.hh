#pragma once
#include "out/LockFreeIO.hh"

namespace labw::art_modern {

class SimpleLFIO : public LockFreeIO<std::string> {
public:
    void write(std::unique_ptr<std::string> value) override { out_ << *value; }
    explicit SimpleLFIO(std::ostream& out)
        : out_(out)
    {
    }

private:
    std::ostream& out_;
};
} // art_modern
// labw
