#pragma once

#define PRINT_MEMORY_USAGE \
        lab::art_modern::details::print_memory_usage(__FILE__, __LINE__);

namespace lab::art_modern::details {
void print_memory_usage(const char* file, int line);


} // lab

