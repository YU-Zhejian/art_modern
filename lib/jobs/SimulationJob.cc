#include "SimulationJob.hh"

namespace labw::art_modern {
SimulationJob::~SimulationJob() { }

SimulationJob::SimulationJob(BaseFastaFetch* fasta_fetch, const CoverageInfo& coverage_info, const int job_id)
    : fasta_fetch(fasta_fetch)
    , coverage_info(coverage_info)
    , job_id(job_id)
{
}
}