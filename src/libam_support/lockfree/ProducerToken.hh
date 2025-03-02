#pragma once

#include <concurrentqueue.h>

namespace labw::art_modern {
    struct ProducerToken{
        moodycamel::ProducerToken token;
    };
} // namespace labw::art_modern
