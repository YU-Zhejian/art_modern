unset(USE_ASIO_PARALLEL)
unset(USE_NOP_PARALLEL)
unset(USE_BS_PARALLEL)

if(NOT DEFINED USE_THREAD_PARALLEL)
    set(USE_THREAD_PARALLEL ASIO)
endif()

if("${USE_THREAD_PARALLEL}" STREQUAL "ASIO")
    set(USE_ASIO_PARALLEL ON)
elseif("${USE_THREAD_PARALLEL}" STREQUAL "NOP")
    set(USE_NOP_PARALLEL ON)
elseif("${USE_THREAD_PARALLEL}" STREQUAL "BS")
    include_directories(BEFORE "${CMAKE_CURRENT_LIST_DIR}/../deps/thread-pool/include/")
    set(USE_BS_PARALLEL ON)
else()
    message(
        FATAL_ERROR
            "Variable USE_THREAD_PARALLEL should be set to one of ASIO, BS or NOP. Current: ${USE_THREAD_PARALLEL}")
endif()
