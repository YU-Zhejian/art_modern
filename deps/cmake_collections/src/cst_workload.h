/*!
 * @file cst_workload.h
 * @author YU Zhejian
 * @brief Workload for testing parallelization libraries like OpenMP, PThread, etc.
 * @version 0.1
 * @date 2024-04-28
 */

#ifndef LIBCONCURRENTQUEUE_CST_WORKLOAD_H
#define LIBCONCURRENTQUEUE_CST_WORKLOAD_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline int perform_sqrt(int num_of_rounds, int num_to_sqrt)
{
    double* result = (double*)malloc(sizeof(double) * num_to_sqrt);
    if (result == NULL) {
        return 1;
    }
    for (int round = 0; round < num_of_rounds; round++) {
        for (int i = 0; i < num_to_sqrt; i++) {
            result[i] = sqrt(i);
        }
    }
    free(result);
    result = NULL;
    return 0;
}

typedef struct
{
    int num_to_sqrt;
    int num_of_rounds;
    int num_of_threads;
} parallel_params_type;

static inline parallel_params_type

parse_args(void)
{
    parallel_params_type parallel_params;
    parallel_params.num_to_sqrt = 10;
    parallel_params.num_of_rounds = 10;
    parallel_params.num_of_threads = 2;
    return parallel_params;
}

#ifdef __cplusplus
};
#endif

#endif // LIBCONCURRENTQUEUE_CST_WORKLOAD_H
