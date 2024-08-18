#pragma once
#include "art_modern_constants.hh"
#include "jobs/SimulationJob.hh"
#include <memory>

namespace labw {
namespace art_modern {

    class JobDispatcher {
    public:
        virtual bool has_next();
        virtual SimulationJob next();
        virtual ~JobDispatcher();
    };

    std::shared_ptr<JobDispatcher> create_job_dispatcher(std::string input_ref_file, std::string input_fcov_param,
        INPUT_FILE_TYPE input_file_type, INPUT_FILE_PARSER input_file_parser, int batch_size);

} // art_modern
} // labw
