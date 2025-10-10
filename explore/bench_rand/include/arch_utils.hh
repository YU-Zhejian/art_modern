#pragma once

#if defined(__ARM_ARCH) || defined(__arm) || defined(__arm__) || defined(__aarch64__)
#define BENCH_RAND_ARCH_ARM
#elif defined(__amd64) || defined(__amd64__) || defined(__x86_64) || defined(__x86_64__)
#define BENCH_RAND_ARCH_X86
#else
#error "Unsupported architecture"
#endif
