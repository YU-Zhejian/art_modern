#pragma once

#include "libam_support/Dtypes.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <string>

namespace labw::art_modern {
class JobExecutor {
public:
    JobExecutor() = default;
    virtual ~JobExecutor() = default;
    DELETE_COPY(JobExecutor)
    DELETE_MOVE(JobExecutor)

    virtual void operator()() = 0;
    [[nodiscard]] virtual bool is_running() const = 0;
    [[nodiscard]] virtual std::string thread_info() const = 0;
};

} // namespace labw::art_modern
