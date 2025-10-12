#pragma once

#include "macros.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <cctype>

namespace Encoder {

    static const std::string base64Chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    static const std::string hexChars =
        "0123456789"
        "ABCDEF";

    static const int8_t hexCharToDec[256] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 0
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 16
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 32
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, -1, -1, -1, -1, -1, // 48
        -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 64
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 80
        -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 96
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 112
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 128
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 144
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 160
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 176
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 192
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 208
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 224
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 240
    };

    static const int8_t b64CharToDec[256] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 0
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 16
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 62, -1, -1, -1, 63, // 32
        52, 53, 54, 55, 56, 57, 58, 59, 60, 61, -1, -1, -1, -1, -1, -1, // 48
        -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, // 64
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, -1, -1, -1, -1, -1, // 80
        -1, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, // 96
        41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, -1, -1, -1, -1, -1, // 112
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 128
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 144
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 160
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 176
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 192
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 208
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 224
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 240
    };

    typedef unsigned char uint8_t;

    inline bool isBase64(uint8_t c) {
        //return (isalnum(c) || (c == '+') || (c == '/'));
        return b64CharToDec[c] != -1;
    }

    inline bool isHex(uint8_t c) {
        //return (isdigit(c) || ((c >= 'a') && (c <= 'f')) || ((c >= 'A') && (c <= 'F')));
        return hexCharToDec[c] != -1;
    }


    inline uint8_t hexToDec(char hex)
    {
        int8_t h = hexCharToDec[(int)hex];
        MYASSERT(h != -1, "not a valid hex character: " << ((int)hex));
        return h;
        //if (hex >= '0' && hex <= '9')
        //    return hex - '0';
        //if (hex >= 'a' && hex <= 'f')
        //    return hex - 'a' + 10;
        //if (hex >= 'A' && hex <= 'F')
        //    return hex - 'A' + 10;
        //char t[2] = { hex, 0 };
        //THROW("not a hex character " << t);
    }

    inline uint8_t b64ToDec(char b64)
    {
        int8_t h = b64CharToDec[(int)b64];
        MYASSERT(h != -1, "not a valid base64 character: " << ((int)b64));
        return h;
        //const size_t nAlpha = 'Z' - 'A' + 1;
        //if (b64 >= 'A') {
        //    if (b64 <= 'Z')
        //        return b64 - 'A';
        //    if (b64 >= 'a' && b64 <= 'z')
        //        return b64 - 'a' + nAlpha;
        //}
        //if (b64 >= '0' && b64 <= '9')
        //    return b64 - '0' + 2 * nAlpha;
        //if (b64 == '+')
        //    return 62;
        //if (b64 == '/')
        //    return 63;
        //char t[2] = { b64, 0 };
        //THROW("not a base 64 character: " << t);
    }

    inline uint8_t hexPairToDec(const char hi, const char lo)
    {
        return (hexToDec(hi) << 4) | hexToDec(lo);
    }


    inline void hexToText(std::string& res, const std::string& hex)
    {
        size_t n = hex.size();
        size_t m = n / 2;
        MYASSERT((n & 0x1) == 0, "n must be even, but got " << n);
        res.resize(m);
        for (size_t i = 0, j = 0; i < m; ++i, j += 2)
            res[i] = hexPairToDec(hex[j], hex[j+1]);
    }

    inline std::ostream& hexToTextStream(std::ostream& res, std::istream& hex)
    {
        char chs[2];
        size_t n = 0;
        for(;;) {
            hex.read(chs, 2);
            n += hex.gcount();
            if (hex.eof())
                break;
            uint8_t txt = hexPairToDec(chs[0], chs[1]);
            res.put(txt);
        }
        MYASSERT(n % 2 == 0, "hex stream length must be even, extracted " << n << " characters");
        return res;
    }

    inline void textToHex(std::string& hex, const std::string& text)
    {
        size_t n = text.size();
        size_t m = n * 2;
        const uint8_t* p = (const uint8_t*) text.c_str();
        hex.resize(m);
        for (size_t i = 0, j = 0; i < n; ++i, j += 2) {
            hex[j] = hexChars[p[i] >> 4];
            hex[j + 1] = hexChars[p[i] & 0xF];
        }
    }


    inline void dec3ToBase64(char res[4], const uint8_t dec[3])
    {
        res[0] = base64Chars[dec[0] >> 2];
        res[1] = base64Chars[((dec[0] & 0x3) << 4) | (dec[1] >> 4)];
        res[2] = base64Chars[((dec[1] & 0xF) << 2) | (dec[2] >> 6)];
        res[3] = base64Chars[dec[2] & 0x3F];
    }

    inline void base64ToDec3(uint8_t res[3], const uint8_t b64[4])
    {
        const uint8_t dec[4] = { b64ToDec(b64[0]), b64ToDec(b64[1]), b64ToDec(b64[2]), b64ToDec(b64[3]) };
        res[0] = (dec[0] << 2) | (dec[1] >> 4);
        res[1] = (dec[1] << 4) | (dec[2] >> 2);
        res[2] = (dec[2] << 6) | dec[3];
    }


    inline void textToBase64(std::string& res, const std::string& txt)
    {
        size_t n = txt.size();
        res.resize((n / 3 + ((n % 3) ? 1 : 0)) * 4);
        const uint8_t* dec = reinterpret_cast<const uint8_t*>(txt.c_str());
        size_t j = 0, i = 0;
        for (; i + 2 < n; i += 3, j += 4)
            dec3ToBase64(&res[j], dec + i);
        switch (n % 3) {
            case 0:
                break;
            case 1:
            {
                uint8_t aux[3] = { dec[i], 0, 0 };
                dec3ToBase64(&res[j], aux);
                res[j + 2] = res[j + 3] = '=';
            }
            break;
            case 2:
            {
                uint8_t aux[3] = { dec[i], dec[i + 1], 0 };
                dec3ToBase64(&res[j], aux);
                res[j + 3] = '=';
            }
            break;
        };
    }

    inline void base64ToText(std::string& res, const std::string& b64)
    {
        size_t n = b64.size();
        MYASSERT(((n % 4) == 0 && n > 0), "bad stream size: " << n);
        size_t nEq = (b64[n - 1] == '=' ? size_t(1) : 0) + (b64[n - 2] == '=' ? size_t(1) : 0);
        size_t n4 = n / 4 - (nEq > 0);
        res.resize(n4 * 3 + ((nEq > 0) ? 3 - nEq : 0));
        const uint8_t* dec = reinterpret_cast<const uint8_t*>(b64.c_str());
        size_t i = 0, j = 0;
        for (; i < n4; ++i, j += 3)
            base64ToDec3(reinterpret_cast<uint8_t*>(&res[j]), dec + 4 * i);
        switch (nEq) {
            case 0:
                break;
            case 1:
            {
                uint8_t r[3];
                const uint8_t b[4] = { dec[4 * i],dec[4 * i + 1],dec[4 * i + 2],'A' };
                base64ToDec3(r, b);
                res[j++] = r[0];
                res[j] = r[1];
            }
            break;
            case 2:
            {
                uint8_t r[3];
                const uint8_t b[4] = { dec[4 * i],dec[4 * i + 1],'A','A' };
                base64ToDec3(r, b);
                res[j] = r[0];
            }
            break;
        };
    }

    inline std::ostream& base64ToTextStream(std::ostream& res, std::istream& b64)
    {
        char b64Ch[4];
        char textCh[3];
        size_t n = 0;
        for (;;) {
            b64.read(b64Ch, 4);
            n += b64.gcount();
            if (b64.eof())
                break;

            if (b64Ch[3] != '=') { // there are no '=' characters
                base64ToDec3((uint8_t*)textCh, (const uint8_t*)b64Ch);
                res.write((char *) textCh, 3);
            }
            else if (b64Ch[2] != '=') { // there is just one '=' character
                b64Ch[3] = 'A';
                base64ToDec3((uint8_t*)textCh, (const uint8_t*)b64Ch);
                res.write(textCh, 2);
            }
            else {  // there are two '=' characters
                b64Ch[2] = b64Ch[3] = 'A';
                base64ToDec3((uint8_t*)textCh, (const uint8_t*)b64Ch);
                res.put(textCh[0]);
            }
        }
        
        MYASSERT(((n % 4) == 0 && n > 0), "base64 stream length must be a multiple of 4, extracted " << n << " characters");
        return res;
    }

};
