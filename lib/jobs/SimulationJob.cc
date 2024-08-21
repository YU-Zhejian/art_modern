#include "SimulationJob.hh"

namespace labw {
namespace art_modern
{
SimulationJob::~SimulationJob()
{
    if(fasta_fetch != nullptr){
        delete fasta_fetch;
    }
}

SimulationJob::SimulationJob(BaseFastaFetch *fasta_fetch, const CoverageInfo &coverage_info, const int job_id)
: fasta_fetch(std::move(fasta_fetch)), coverage_info(coverage_info), job_id(job_id)
{

}
}}