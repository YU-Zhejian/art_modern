#include <libdeflate.h>

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    printf("libdeflate version: %s\n", LIBDEFLATE_VERSION_STRING);
    // Generate a random string of 4k chars
    const size_t data_size = 4096;
    char *data = (char *)malloc(data_size);
    if (data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }
    for (size_t i = 0; i < data_size - 1; i++) {
        data[i] = 'A' + (rand() % 26);
    }
    data[data_size - 1] = '\0';
    // Compress the data
    struct libdeflate_compressor *compressor = libdeflate_alloc_compressor(6);
    if (compressor == NULL) {
        fprintf(stderr, "Failed to allocate compressor\n");
        free(data);
        return EXIT_FAILURE;
    }
    size_t max_compressed_size = libdeflate_deflate_compress_bound(compressor, data_size);
    void *compressed_data = malloc(max_compressed_size);
    if (compressed_data == NULL) {
        fprintf(stderr, "Memory allocation for compressed data failed\n");
        libdeflate_free_compressor(compressor);
        free(data);
        return EXIT_FAILURE;
    }
    size_t compressed_size = libdeflate_deflate_compress(compressor, data, data_size, compressed_data, max_compressed_size);
    libdeflate_free_compressor(compressor);
    if (compressed_size == 0) {
        fprintf(stderr, "Compression failed\n");
        free(compressed_data);
        free(data);
        return EXIT_FAILURE;
    }
    // Decompress the data
    struct libdeflate_decompressor *decompressor = libdeflate_alloc_decompressor();
    if (decompressor == NULL) {
        fprintf(stderr, "Failed to allocate decompressor\n");
        free(compressed_data);
        free(data);
        return EXIT_FAILURE;
    }
    void *decompressed_data = malloc(data_size);
    if (decompressed_data == NULL) {
        fprintf(stderr, "Memory allocation for decompressed data failed\n");
        libdeflate_free_decompressor(decompressor);
        free(compressed_data);
        free(data);
        return EXIT_FAILURE;
    }
    enum libdeflate_result result = libdeflate_deflate_decompress(decompressor, compressed_data, compressed_size, decompressed_data, data_size, NULL);
    libdeflate_free_decompressor(decompressor);
    if (result != LIBDEFLATE_SUCCESS) {
        fprintf(stderr, "Decompression failed: %d\n", result);
        free(decompressed_data);
        free(compressed_data);
        free(data);
        return EXIT_FAILURE;
    }
    // Verify the decompressed data matches the original data
    if (memcmp(data, decompressed_data, data_size) != 0) {
        fprintf(stderr, "Decompressed data does not match original data\n");
        free(decompressed_data);
        free(compressed_data);
        free(data);
        return EXIT_FAILURE;
    }
    printf("Compression and decompression successful!\n");
    free(decompressed_data);
    free(compressed_data);
    free(data);
    return EXIT_SUCCESS;
}
