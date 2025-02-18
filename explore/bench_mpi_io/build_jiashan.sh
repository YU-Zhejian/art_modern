#!/usr/bin/env bash
module () 
{ 
    eval `/usr/bin/modulecmd bash $*`
}
module load mpi/openmpi/4.1.4-gcc

make build CMAKE_FLAGS="-DCMAKE_CXX_COMPILER=$(which g++)"
