if(NOT DEFINED CEU_CM_GIT_ADAPTER_WAS_ALREADY_INCLUDED)
    set(CEU_CM_GIT_ADAPTER_WAS_ALREADY_INCLUDED
        TRUE
        CACHE INTERNAL "Whether the CEU CM Git adapter was already included.")
    set(CEU_CM_GIT_COMMIT_HASH
        "Do not have Git VCS"
        CACHE STRING "Default to unfound")
    set(CEU_CM_GIT_COMMIT_DATE
        "Do not have Git VCS"
        CACHE STRING "Default to unfound")
    set(CEU_CM_GIT_GIT_COMMIT_MESSAGE
        "Do not have Git VCS"
        CACHE STRING "Default to unfound")

    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" log -1 --pretty=format:%H
            OUTPUT_VARIABLE CEU_CM_GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" log -1 --pretty=format:%aD
            OUTPUT_VARIABLE CEU_CM_GIT_COMMIT_DATE
            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" log -1 --pretty=format:%s
            OUTPUT_VARIABLE CEU_CM_GIT_GIT_COMMIT_MESSAGE
            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    else()
        message(WARNING "Git executable not found.")
    endif()
    message(STATUS "/------------------- Finding and assessing git -------------------\\")
    message(STATUS "|git executable found at ${GIT_EXECUTABLE}")
    message(
        STATUS
            "|On Git commit ${CEU_CM_GIT_GIT_COMMIT_MESSAGE} (${CEU_CM_GIT_COMMIT_HASH}) at ${CEU_CM_GIT_COMMIT_DATE}")
    message(STATUS "\\------------------- Finding and assessing git -------------------/")
endif()
