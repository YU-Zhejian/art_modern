unset(USE_ASIO_PARALLEL)
unset(USE_NOP_PARALLEL)

if("${USE_THREAD_PARALLEL}" STREQUAL "ASIO")
    set(USE_ASIO_PARALLEL ON)
elseif("${USE_THREAD_PARALLEL}" STREQUAL "NOP")
    set(USE_NOP_PARALLEL ON)
else()
    message(
        FATAL_ERROR "Variable USE_THREAD_PARALLEL should be set to one of ASIO or NOP. Current: ${USE_THREAD_PARALLEL}")
endif()
