#pragma once
#include <string>
#include "BaseReadOutput.hh"

namespace labw::art_modern {


class BaseFileReadOutput  : public BaseReadOutput {
public:
    explicit BaseFileReadOutput(const std::string& filename);
    void close() override;

protected:
    const std::string filename;
    bool is_closed_;
};

}
