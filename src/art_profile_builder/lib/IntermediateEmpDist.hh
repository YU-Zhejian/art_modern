#pragma once

#include "art_profile_builder/lib/IntermediateEmpDistPosition.hh"

#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/sam.h>

#include <cstdlib>
#include <iostream>
#include <vector>

namespace labw::art_modern {
class IntermediateEmpDist {
public:
    explicit IntermediateEmpDist(std::size_t read_length);
    ~IntermediateEmpDist() = default;

    IntermediateEmpDist(const IntermediateEmpDist& ied) = default;
    IntermediateEmpDist& operator=(const IntermediateEmpDist&) = delete;

    DELETE_MOVE(IntermediateEmpDist)

    bool parse_read(const bam1_t* b);

    void accumulate();

    void add(const IntermediateEmpDist& other);

    void write(std::ostream& oss, bool is_ob) const;

private:
    std::vector<IntermediateEmpDistPosition> positions_;
    const std::size_t read_length_;
};
} // namespace labw::art_modern
