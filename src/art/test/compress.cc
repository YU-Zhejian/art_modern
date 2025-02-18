// Implement LZW
// Implement BWT+RLE

#include "art/builtin_profiles.hh"

#include <pcg_random.hpp>

#include <cstddef>
#include <stdexcept>
#include <vector>
#include <string>
#include <array>
#include <random>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <cstdint>

using namespace labw::art_modern;

class BWTFormat{
    constexpr static const char TERM = 0x7F;
public:
    static std::string compress(const std::string& src){
        // Scan string for TERM and except if occurs
        for (char const c : src) {
            if (c == TERM) {
                throw std::runtime_error("Input string contains TERM character");
            }
        }
        std::string const s = src + TERM; // Append the delimiter
        std::vector<std::string> rotations;

        // Generate all rotations of the string
        for (size_t i = 0; i < s.length(); ++i) {
            rotations.push_back(s.substr(i) + s.substr(0, i));
        }

        // Sort the rotations lexicographically
        std::sort(rotations.begin(), rotations.end());

        // Extract the last column
        std::string bwt;
        for (const auto& rotation : rotations) {
            bwt += rotation.back();
        }

        return bwt;
    }
    static std::string decompress(const std::string& src){
        std::vector<std::string> table(src.length());

        // Reconstruct the table
        for (size_t i = 0; i < src.length(); ++i) {
            for (size_t j = 0; j < src.length(); ++j) {
                table[j] = src[j] + table[j];
            }
            std::sort(table.begin(), table.end());
        }

        // Find the row that ends with the delimiter
        for (const auto& row : table) {
            if (row.back() == TERM) {
                return row.substr(0, row.length() - 1); // Remove the delimiter
            }
        }

        return "";
    }
};

class RLEFormat{
private:
    std::vector<uint32_t> runLengths;
    std::vector<char> runValues;
    RLEFormat() = default;
    RLEFormat(const std::vector<uint32_t>& runLengths, const std::vector<char>& runValues) : runLengths(runLengths), runValues(runValues) {}
    explicit RLEFormat(const std::string& src){
        char currentChar = src[0];
        uint32_t count = 1;

        for (size_t i = 1; i < src.size(); ++i) {
            if (src[i] == currentChar) {
                count++;
            } else {
                runValues.push_back(currentChar);
                runLengths.push_back(count);  // Store count as char
                currentChar = src[i];
                count = 1;
            }
        }

        // Add the last character and its count
        runValues.push_back(currentChar);
        runLengths.push_back(count);
    }
    std::string encode(){
        std::ostringstream ss;
        uint32_t len = runLengths.size();
        ss << std::string (reinterpret_cast<const char*>(&len), 4);
        for (auto rl : runLengths){
            ss << std::string (reinterpret_cast<const char*>(&rl), 4);
        }
        for (auto rv : runValues){
            ss << rv;
        }
        return ss.str();
    }
    static RLEFormat decode(const std::string& encoded){
        uint32_t const len = *reinterpret_cast<const uint32_t*>(encoded.data());

        std::vector<uint32_t> runLengths;
        for (size_t i = 0; i < len; ++i) {
            runLengths.push_back(*reinterpret_cast<const uint32_t*>(encoded.data()+4+i*4));
        }
        std::vector<char> const runValues = std::vector<char>(encoded.data()+4+len*4, encoded.data()+encoded.size());
        return RLEFormat{runLengths, runValues};
    }
    std::string to_original(){
        std::ostringstream ss;
        for (size_t i = 0; i < runLengths.size(); ++i) {
            for (uint32_t j = 0; j < runLengths[i]; ++j) {
                ss << runValues[i];
            }
        }
        return ss.str();
    }
public:
    static std::string compress(const std::string& src){
        RLEFormat rle = RLEFormat(src);
        return rle.encode();
    }
    static std::string decompress(const std::string& src){
        RLEFormat rle = RLEFormat::decode(src);
        return rle.to_original();
    }
};

int main() {
    std::string const input = "Hello, World!";
    std::cout << "RLE C->DC: " << RLEFormat::decompress(RLEFormat::compress(input)) << std::endl;
    std::cout << "BWT C->DC: " << BWTFormat::decompress(BWTFormat::compress(input)) << std::endl;
    pcg32_fast rng{std::random_device{}()};
    std::uniform_int_distribution<int>dist{0, N_BUILTIN_PROFILE - 1};
    for (int i = 0; i < 1000; i++) {
        const std::string str{ENCODED_BUILTIN_PROFILES[dist(rng)][0]};
        std::cout << BWTFormat::compress(str);
        const std::string bwt_rle_cmp = RLEFormat::compress(BWTFormat::compress(str));
        const std::string bwt_rle_restored = BWTFormat::decompress(RLEFormat::decompress(bwt_rle_cmp));
        if(str != bwt_rle_restored){
            std::cout << "Error: " << str << " != " << bwt_rle_restored << std::endl;
            return 1;
        }
        const double cmp_ratio = static_cast<double>(bwt_rle_cmp.size()) / str.size();
        std::cout << "Compression ratio: " << cmp_ratio << std::endl;
    }
    return 0;
}
