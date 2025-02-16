#include "art_modern_capi.h"
#include "art_modern_config.h"

#include "art/lib/ArtConstants.hh"
#include "art/lib/ArtJobExecutor.hh"
#include "art/lib/ArtParams.hh"

#include "libam_support/Constants.hh"
#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/jobs/JobPool.hh"
#include "libam_support/jobs/SimulationJob.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutputDispatcher.hh"
#include "libam_support/ref/batcher/FastaStreamBatcher.hh"
#include "libam_support/ref/batcher/InMemoryFastaBatcher.hh"
#include "libam_support/ref/batcher/Pbsim3TranscriptBatcher.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/ref/fetch/FaidxFetch.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"

#include <boost/log/trivial.hpp>

#include <atomic>
#include <chrono>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>
#include <utility>

using namespace labw::art_modern;

amcapi_params_t* amcapi_init_params(void)
{
    auto* retv = new amcapi_params_t();
    retv->am_params = new ArtParams {};
    return retv;
}

amcapi_pwa_t* amcapi_generate(void) { }