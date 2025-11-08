if(WITH_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_C_HEADER_DIR} ${MPI_CXX_HEADER_DIR})
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} MPI::MPI_C) # No need to add MPI::MPI_CXX
    set(WITH_MPI ON)
endif()
