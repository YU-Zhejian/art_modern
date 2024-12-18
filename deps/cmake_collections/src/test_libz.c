#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

uLong get_adler_crc32(char* buff)
{
    uLong adler_checksum = adler32(0L, Z_NULL, 0);
    adler_checksum = adler32(adler_checksum, buff, sizeof(buff));
    return adler_checksum;
}

int main(void)
{
    printf("Zlib version: %s\n", ZLIB_VERSION);

    char strSrc[] = "hello world!";
    char buf[1024] = { 0 };
    char strDst[1024] = { 0 };
    uLong srcLen = sizeof(strSrc);
    uLong bufLen = sizeof(buf);
    uLong dstLen = sizeof(strDst);

    printf("Src string: %s (len = %zu)\n", strSrc, srcLen);

    compress(buf, &bufLen, strSrc, srcLen);
    printf("After Compressed Length: %zu\n", bufLen);

    uncompress(strDst, &dstLen, buf, bufLen);
    printf("Dst string: %s (len = %zu)\n", strDst, dstLen);

    if (strncmp(strSrc, strDst, 1024) != 0) {
        return 1;
    }
    if (get_adler_crc32(strSrc) != get_adler_crc32(strDst)) {
        return 1;
    }
    return 0;
}
