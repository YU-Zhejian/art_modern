/*

  https://github.com/superwills/NibbleAndAHalf
  base64.h -- Fast base64 encoding and decoding.
  version 1.0.0, April 17, 2013 143a

  Copyright (C) 2013 William Sherif

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  William Sherif
  will.sherif@gmail.com

  YWxsIHlvdXIgYmFzZSBhcmUgYmVsb25nIHRvIHVz

  YU Zhejian fixed some clang-tidy warnings.
*/
#ifndef BASE64_H
#define BASE64_H

#include <cstddef>
#include <string>

constexpr static char b64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// maps A=>0,B=>1..
constexpr static unsigned char unb64[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //10
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //20
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //30
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //40
        0, 0, 0, 62, 0, 0, 0, 63, 52, 53, //50
        54, 55, 56, 57, 58, 59, 60, 61, 0, 0, //60
        0, 0, 0, 0, 0, 0, 1, 2, 3, 4, //70
        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, //80
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, //90
        25, 0, 0, 0, 0, 0, 0, 26, 27, 28, //100
        29, 30, 31, 32, 33, 34, 35, 36, 37, 38, //110
        39, 40, 41, 42, 43, 44, 45, 46, 47, 48, //120
        49, 50, 51, 0, 0, 0, 0, 0, 0, 0, //130
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //140
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //150
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //160
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //170
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //180
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //190
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //200
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //210
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //220
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //230
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //240
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //250
        0, 0, 0, 0, 0, 0,
}; // This array has 256 elements

// Converts binary data of length=len to base64 characters.
[[maybe_unused]] static inline std::string base64_encode(const void *data, std::size_t length) {
    const auto *bin = static_cast<const unsigned char *>(data);

    auto modLength = length % 3;
    // 2 gives 1 and 1 gives 2, but 0 gives 0.
    auto padding = ((modLength & 1U) << 1U) + ((modLength & 2U) >> 1U);

    std::string res;
    res.reserve(4 * (length + padding) / 3);

    std::size_t byteNo = 0;
    for (; byteNo <= length - 3; byteNo += 3) {
        const unsigned char BYTE0 = bin[byteNo];
        const unsigned char BYTE1 = bin[byteNo + 1];
        const unsigned char BYTE2 = bin[byteNo + 2];
        res.append(1, b64[BYTE0 >> 2U]);
        res.append(1, b64[((0x3 & BYTE0) << 4U) + (BYTE1 >> 4U)]);
        res.append(1, b64[((0x0f & BYTE1) << 2U) + (BYTE2 >> 6U)]);
        res.append(1, b64[0x3f & BYTE2]);
    }

    if (padding == 2) {
        res.append(1, b64[bin[byteNo] >> 2U]);
        res.append(1, b64[(0x3 & bin[byteNo]) << 4U]);
        res.append(1, '=');
        res.append(1, '=');
    } else if (padding == 1) {
        res.append(1, b64[bin[byteNo] >> 2U]);
        res.append(1, b64[((0x3 & bin[byteNo]) << 4) + (bin[byteNo + 1] >> 4U)]);
        res.append(1, b64[(0x0f & bin[byteNo + 1]) << 2]);
        res.append(1, '=');
    }
    return res;
}

static inline std::string base64_decode(const char *base64, std::size_t length) {
    const auto *data = reinterpret_cast<const unsigned char *>(base64);
    // 2 accesses below would be OOB.
    if (length < 2) {
        return "";
    }

    int padding = 0;
    if (data[length - 1] == '=') { ++padding; }
    if (data[length - 2] == '=') { ++padding; }

    std::string res;
    res.reserve(3 * (length >> 2U) - padding);

    std::size_t charNo = 0;
    for (; charNo <= length - 4 - padding; charNo += 4) {
        const unsigned char A = unb64[data[charNo]];
        const unsigned char B = unb64[data[charNo + 1]];
        const unsigned char C = unb64[data[charNo + 2]];
        const unsigned char D = unb64[data[charNo + 3]];

        res.append(1, static_cast<char>((A << 2U) | (B >> 4U)));
        res.append(1, static_cast<char>((B << 4U) | (C >> 2U)));
        res.append(1, static_cast<char>((C << 6U) | (D)));
    }

    if (padding == 1) {
        const unsigned char A = unb64[data[charNo]];
        const unsigned char B = unb64[data[charNo + 1]];
        const unsigned char C = unb64[data[charNo + 2]];

        res.append(1, static_cast<char>((A << 2U) | (B >> 4U)));
        res.append(1, static_cast<char>((B << 4U) | (C >> 2U)));
    } else if (padding == 2) {
        const unsigned char A = unb64[data[charNo]];
        const unsigned char B = unb64[data[charNo + 1]];

        res.append(1, static_cast<char>((A << 2U) | (B >> 4U)));
    }

    return res;
}

#endif
