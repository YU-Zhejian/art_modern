#include "cst_workload.h"

#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

void worker(int id, int num_to_sqrt, int num_of_rounds, int* retv)
{
    printf("Thread %d start\n", id);
    *retv = perform_sqrt(num_of_rounds, num_to_sqrt);
    printf("Thread %d join with return value %d\n", id, *retv);
}

int main()
{
    // Print ID of this thread
    const auto hash_this_thread = std::hash<std::thread::id> {}(std::this_thread::get_id());
    printf("Main thread ID: %zu\n", hash_this_thread);
    parallel_params_type parallel_params = parse_args();
    std::vector<std::thread> thread_vector;
    thread_vector.reserve(parallel_params.num_of_threads);
    std::vector<int> return_vector(parallel_params.num_of_threads);
    for (int i = 0; i < parallel_params.num_of_threads; i++) {
        thread_vector.emplace_back(worker, i, parallel_params.num_to_sqrt, parallel_params.num_of_rounds,
            &return_vector[i]);
    }
    for (int i = 0; i < parallel_params.num_of_threads; i++) {
        if (!thread_vector[i].joinable()) {
            printf("Thread %d not joinable\n", i);
            return EXIT_FAILURE;
        }
        thread_vector[i].join();
        printf("Thread %d join with return value %d\n", i, return_vector[i]);
    }
    return EXIT_SUCCESS;
}
