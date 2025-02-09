#pragma once

// NOLINTBEGIN
#define PRINT_MEMORY_USAGE lab::art_modern::details::print_memory_usage(__FILE__, __LINE__);
// NOLINTEND

namespace lab::art_modern::details {
[[maybe_unused]] void print_memory_usage(const char* file, int line);

} // namespace lab::art_modern::details
