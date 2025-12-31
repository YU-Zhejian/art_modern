#pragma once

#include <cstdint>
namespace labw::art_modern {
enum class APB_FORMAT : std::uint8_t { FASTQ, SAM, BAM, CRAM, SRA, AUTO };
constexpr char APB_FORMAT_FASTQ_STR[] = "FASTQ";
constexpr char APB_FORMAT_SAM_STR[] = "SAM";
constexpr char APB_FORMAT_BAM_STR[] = "BAM";
constexpr char APB_FORMAT_CRAM_STR[] = "CRAM";
constexpr char APB_FORMAT_SRA_STR[] = "SRA";
constexpr char APB_FORMAT_AUTO_STR[] = "AUTO";
constexpr char SRA_SIGNATURE[8] = { 'N', 'C', 'B', 'I', '.', 's', 'r', 'a' };
} // namespace labw::art_modern
