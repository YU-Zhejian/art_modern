#ifndef _TWOBIT_H
#define _TWOBIT_H

#include <stdint.h>
#include <stdlib.h>

typedef struct twobit_decoder {
    char* BYTE_TO_BASE_PRE_COMPUTED;
} twobit_decoder_t;

twobit_decoder_t* twobit_decoder_init(void);
void twobit_decoder_free(twobit_decoder_t* twobit_decoder);

void twobit_to_nt(twobit_decoder_t* twobit_decoder, const uint8_t* const uint8s, char* nts, size_t src_pos,
    size_t dst_pos, size_t num_bytes_to_read);
#endif
