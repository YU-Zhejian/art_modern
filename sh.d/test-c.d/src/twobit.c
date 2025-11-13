#include "twobit.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

static const char* const BYTE_TO_BASE = "TCAG";

twobit_decoder_t* twobit_decoder_init(void)
{
    twobit_decoder_t* twobit_decoder = (twobit_decoder_t*)malloc(sizeof(twobit_decoder_t));
    twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED = (char*)calloc((unsigned long)(256 * 4), sizeof(char));
    for (int i = 0; i < 256; i++) {
        twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED[(i << 2) + 0] = BYTE_TO_BASE[(i >> 6) & 3];
        twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED[(i << 2) + 1] = BYTE_TO_BASE[(i >> 4) & 3];
        twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED[(i << 2) + 2] = BYTE_TO_BASE[(i >> 2) & 3];
        twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED[(i << 2) + 3] = BYTE_TO_BASE[i & 3];
    }
    return twobit_decoder;
}

void twobit_decoder_free(twobit_decoder_t* twobit_decoder)
{
    free(twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED);
    free(twobit_decoder);
}

void twobit_to_nt(twobit_decoder_t* twobit_decoder, const uint8_t* const uint8s, char* nts, size_t src_pos, size_t dst_pos,
    size_t num_bytes_to_read)
{
    for (size_t i = 0; i < num_bytes_to_read; i++) {
        uint8_t byte = uint8s[src_pos++];
        nts[dst_pos++] = twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED[(byte << 2) + 0];
        nts[dst_pos++] = twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED[(byte << 2) + 1];
        nts[dst_pos++] = twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED[(byte << 2) + 2];
        nts[dst_pos++] = twobit_decoder->BYTE_TO_BASE_PRE_COMPUTED[(byte << 2) + 3];
    }
}
