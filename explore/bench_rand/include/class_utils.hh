#pragma once

// Define a macro that deletes copy and move constructors and assignment operators
#define DELETE_COPY_MOVE(ClassName)                                                                                    \
    ClassName(const ClassName&) = delete;                                                                              \
    ClassName(ClassName&&) = delete;                                                                                   \
    ClassName& operator=(const ClassName&) = delete;                                                                   \
    ClassName& operator=(ClassName&&) = delete;

#define SPAN_RESULT_TYPE                                                                                               \
    static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }                             \
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }