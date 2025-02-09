#ifndef AMCAPI_H
#define AMCAPI_H

#include <stdbool.h>
#include <stdint.h>

#include <htslib/hts.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /**
     * Query sequence with gap inserted using -
     */
    char* aln_query;
    /**
     * Reference sequence with gap inserted using -
     */
    char* aln_ref;
    /**
     * Query sequence without gap.
     */
    char* query;
    /**
     * Reference sequence without gap.
     */
    char* ref;
    /**
     * Quality sequence whose elngth should be the same as query.
     */
    char* qual;
    char* read_name;
    char* contig_name;
    hts_pos_t pos_on_contig;
    bool is_plus_strand;
} amcapi_pwa_t;

typedef struct {
    void* am_params;
} amcapi_params_t;

amcapi_params_t* amcapi_init_params(void);

amcapi_pwa_t* amcapi_generate(void);

#ifdef __cplusplus
}
#endif
#endif
