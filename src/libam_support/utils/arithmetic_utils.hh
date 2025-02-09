#pragma once

namespace labw::art_modern {
template <typename T> T am_max(T a, T b)
{
    if (a > b) {
        return a;
    }
    return b;
}
template <typename T> T am_min(T a, T b)
{
    if (a < b) {
        return a;
    }
    return b;
}
} // namespace labw::art_modern