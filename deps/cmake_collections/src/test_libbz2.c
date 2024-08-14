//
// Created by yuzj on 3/23/22.
//
#include <bzlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    printf("libbz2 version: %s\n", BZ2_bzlibVersion());

    char strSrc[] = "hello world!";
    char buf[1024] = { 0 };
    char strDst[1024] = { 0 };
    unsigned srcLen = sizeof(strSrc);
    unsigned bufLen = sizeof(buf);
    unsigned dstLen = sizeof(strDst);

    printf("Src string: %s (len = %u)\n", strSrc, srcLen);

    BZ2_bzBuffToBuffCompress(buf, &bufLen, strSrc, srcLen, 1, 0, 0);
    printf("After Compressed Length: %u\n", bufLen);

    BZ2_bzBuffToBuffDecompress(strDst, &dstLen, buf, bufLen, 0, 0);
    printf("Dst string: %s (len = %u)\n", strDst, dstLen);

    if (strcmp(strSrc, strDst) != 0) {
        return 1;
    }
    return 0;
}
