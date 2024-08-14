#include "cst_workload.h"
#include <pthread.h>
#include <stdio.h>

typedef struct
{
    int thread_id;
    int num_to_sqrt;
    int num_of_rounds;
    int retv;
} pthread_params_type;

void* pthread_thread(void* args)
{
    pthread_params_type* pthread_params = (pthread_params_type*)args;
    printf("Thread %d start\n", pthread_params->thread_id);
    pthread_params->retv = perform_sqrt(pthread_params->num_of_rounds, pthread_params->num_to_sqrt);
    printf("Thread %d end with return value %d\n", pthread_params->thread_id, pthread_params->retv);
    pthread_exit(NULL);
}

int main(void)
{
    parallel_params_type parallel_params = parse_args();
    pthread_t* thr = (pthread_t*)malloc(sizeof(pthread_t) * parallel_params.num_of_threads);
    pthread_params_type* pthread_params_arr = (pthread_params_type*)malloc(sizeof(pthread_params_type) * parallel_params.num_of_threads);
    for (int i = 0; i < parallel_params.num_of_threads; ++i) {
        pthread_params_arr[i].thread_id = i;
        pthread_params_arr[i].num_to_sqrt = parallel_params.num_to_sqrt;
        pthread_params_arr[i].num_of_rounds = parallel_params.num_of_rounds;
        pthread_create(&thr[i], NULL, pthread_thread, (void*)&pthread_params_arr[i]);
    }
    for (int i = 0; i < parallel_params.num_of_threads; ++i) {
        pthread_join(thr[i], NULL);
        printf("Thread %d join with return value %d\n", i, pthread_params_arr[i].retv);
    }
    free(pthread_params_arr);
    free(thr);
    return 0;
}
