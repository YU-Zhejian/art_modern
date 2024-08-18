#include "JobDispatcher.hh"
#include "NotImplementedException.hh"

namespace labw {
namespace art_modern {
    bool JobDispatcher::has_next() { throw NotImplementedException(); }
    SimulationJob JobDispatcher::next() { throw NotImplementedException(); }
    JobDispatcher::~JobDispatcher() = default;

} // art_modern
} // labw