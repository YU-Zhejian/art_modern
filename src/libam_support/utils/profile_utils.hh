/**
 * @file profile_utils.hh
 * @brief Utility functions for profiling and memory usage tracking.
 *
 * This file contains macros and function declarations to assist in profiling
 * and tracking memory usage in the application. The PRINT_MEMORY_USAGE macro
 * can be used to log memory usage at specific points in the code.
 */
#pragma once

// NOLINTBEGIN
#define PRINT_MEMORY_USAGE lab::art_modern::details::print_memory_usage(__FILE__, __LINE__);
// NOLINTEND

namespace lab::art_modern::details {
[[maybe_unused]] void print_memory_usage(const char* file, int line);

} // namespace lab::art_modern::details
