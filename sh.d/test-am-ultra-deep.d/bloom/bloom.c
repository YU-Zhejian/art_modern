#include "bloom.h"
#include "../dtypes.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct bloom_hash {
    hash_function func;
    struct bloom_hash* next;
};

struct bloom_filter {
    struct bloom_hash* func;
    void* bits;
    size_t size;
};

bloom_t bloom_create(size_t size)
{
    bloom_t res = calloc(1, sizeof(struct bloom_filter));
    res->size = size;
    res->bits = malloc(size);
    memset(res->bits, 0, size);
    return res;
}

void bloom_free(bloom_t filter)
{
    if (filter) {
        while (filter->func) {
            struct bloom_hash* h = filter->func;
            filter->func = h->next;
            free(h);
        }
        free(filter->bits);
        free(filter);
    }
}

void bloom_add_hash(bloom_t filter, hash_function func)
{
    struct bloom_hash* h = calloc(1, sizeof(struct bloom_hash));
    h->func = func;
    struct bloom_hash* last = filter->func;
    while (last && last->next) {
        last = last->next;
    }
    if (last) {
        last->next = h;
    } else {
        filter->func = h;
    }
}

void bloom_add(bloom_t filter, const void* item, const size_t len)
{
    struct bloom_hash* h = filter->func;
    uint8_t* bits = filter->bits;
    while (h) {
        hash_type hash = h->func(item, len);
        hash %= filter->size * 8;
        bits[hash / 8] |= 1 << hash % 8;
        h = h->next;
    }
}

bool bloom_test(bloom_t filter, const void* item, const size_t len)
{
    struct bloom_hash* h = filter->func;
    uint8_t* bits = filter->bits;
    while (h) {
        hash_type hash = h->func(item, len);
        hash %= filter->size * 8;
        if (!(bits[hash / 8] & 1 << hash % 8)) {
            return false;
        }
        h = h->next;
    }
    return true;
}
