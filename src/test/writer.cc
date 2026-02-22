/**
 * Copyright 2026 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#define BOOST_TEST_MODULE test_writer // NOLINT

#include "art_modern_config.h" // NOLINT: For WITH_LIBDEFLATE

#include <boost/filesystem/operations.hpp>

#include "libam_support/writer/BGZipWriter.hh"
#include "libam_support/writer/SimpleWriter.hh"
#include "libam_support/writer/WriterInterface.hh"
#if 0
#include "libam_support/writer/ZlibWriter.hh"
#ifdef WITH_LIBDEFLATE
#include "libam_support/writer/LibDeflateWriter.hh"
#endif
#endif

#include "libam_support/ds/pcg_32_c.hh"

#include <boost/log/trivial.hpp>
#include <boost/test/unit_test.hpp>

#include <zlib.h>

#include <numeric>
#include <random>
#include <string>

using namespace labw::art_modern;

namespace {
std::string read_gzip(const std::string& file_name)
{
    gzFile gz_file = gzopen(file_name.c_str(), "rb");
    if (gz_file == nullptr) {
        throw std::runtime_error("Failed to open gzip file: " + file_name);
    }
    std::string result;
    char buffer[4096];
    int bytes_read = 0;
    while ((bytes_read = gzread(gz_file, buffer, sizeof(buffer))) > 0) {
        result.append(buffer, bytes_read);
    }
    gzclose(gz_file);
    return result;
}

void assess(const std::string& data, std::unique_ptr<WriterInterface> writer,
    const std::function<std::string(const std::string&)>& read_fn)
{
    const auto file_name = writer->get_filename();
    writer->write(data);
    writer->flush();
    writer->close();
    // Use Boost to tell actual file size
    const auto read_data = read_fn(file_name);
    BOOST_TEST(read_data == data, "" << writer->name() << " read data does not match written data!");
}
} // namespace

BOOST_AUTO_TEST_CASE(test_writer)
{
    // Make temp dir
    const auto temp_dir = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path();
    boost::filesystem::create_directory(temp_dir);
    // Generate 2K data
    pcg32_c rng(42);
    std::uniform_int_distribution<char> dist(std::numeric_limits<char>::min(), std::numeric_limits<char>::max());
    std::string data(2048, '\0');
    for (auto& c : data) {
        c = dist(rng);
    }
    for (const auto& buffer_size : { 0, 1, 1024, 2048, 4096 }) {
        assess(data,
            std::make_unique<SimpleWriter>(
                (temp_dir / ("simple_" + std::to_string(buffer_size) + ".txt")).string(), buffer_size),
            [](const std::string& path) {
                std::ifstream in(path);
                return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            });
        for (int compression_level = 0; compression_level <= 9; ++compression_level) {
#if 0
            assess(data, std::make_unique<ZLibWriter>((temp_dir / ("zlib_" + std::to_string(buffer_size) + ".gz")).string(), compression_level, buffer_size), read_gzip);
#endif
            for (auto const num_threads : { 1, 2, 4 }) {
                assess(data,
                    std::make_unique<BGZipWriter>(
                        (temp_dir
                            / ("bgzip_" + std::to_string(buffer_size) + "_" + std::to_string(compression_level) + "_"
                                + std::to_string(num_threads) + ".bgz"))
                            .string(),
                        compression_level, buffer_size, false, num_threads),
                    read_gzip);
                assess(data,
                    std::make_unique<BGZipWriter>(
                        (temp_dir
                            / ("bgzip_" + std::to_string(buffer_size) + "_" + std::to_string(compression_level) + "_"
                                + std::to_string(num_threads) + ".gz"))
                            .string(),
                        compression_level, buffer_size, true, num_threads),
                    read_gzip);
            }
        }
#if 0
#ifdef WITH_LIBDEFLATE
        for (int compression_level = 0; compression_level <= 12; ++compression_level)
        {
            assess(data, std::make_unique<LibDeflateWriter>((temp_dir / ("libdeflate_" + std::to_string(buffer_size) + "_" + std::to_string(compression_level) + ".deflate")).string(), compression_level, buffer_size), read_gzip);
        }
#endif
#endif
        // Delete temp dir
    }
    boost::filesystem::remove_all(temp_dir);
}
