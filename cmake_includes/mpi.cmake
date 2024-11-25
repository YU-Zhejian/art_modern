if(WITH_MPI)
    find_package(MPI REQUIRED) # TODO
    include_directories(${MPI_C_HEADER_DIR} ${MPI_CXX_HEADER_DIR})
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} MPI::MPI_C MPI::MPI_CXX)
    set(WITH_MPI ON)
    set(WITH_PROTOBUF ON)
endif()
