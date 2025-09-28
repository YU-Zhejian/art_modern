#pragma once
#include "class_utils.hh"

#include <absl/base/attributes.h>
#include <random>

template <typename xoroshiro_type, typename result_type, int num_states = 2> class XoroshiroWrapper {
public:
    XoroshiroWrapper()
    {
        std::random_device rd;
        for (int i = 0; i < num_states; ++i) {
            xoroshiro_impl_.s[i] = rd();
        }
    }
    DELETE_COPY_MOVE(XoroshiroWrapper)
    ~XoroshiroWrapper() = default;
    ABSL_ATTRIBUTE_ALWAYS_INLINE result_type operator()() { return xoroshiro_impl_.next(); }
    SPAN_RESULT_TYPE

private:
    xoroshiro_type xoroshiro_impl_ {};
};
