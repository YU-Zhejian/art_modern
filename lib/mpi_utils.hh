#pragma once
#include <cstdlib>
namespace labw {
namespace art_modern {

    void exit_mpi(const int status);
    void abort_mpi(const int status = EXIT_FAILURE);
}
}
