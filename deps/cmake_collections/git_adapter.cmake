# FIXME: Code not working.
if(NOT DEFINED CEU_CM_GIT_ADAPTER_WAS_ALREADY_INCLUDED)
    set(CEU_CM_GIT_ADAPTER_WAS_ALREADY_INCLUDED
        ON
        CACHE INTERNAL "Whether the CEU CM Git adapter was already included.")
    set(CEU_CM_GIT_COMMIT_HASH
        "UNKNOWN"
        CACHE STRING "Default if Git was not found")
    set(CEU_CM_GIT_COMMIT_DATE
        "UNKNOWN"
        CACHE STRING "Default if Git was not found")
    set(CEU_CM_GIT_COMMIT_MESSAGE
        "UNKNOWN"
        CACHE STRING "Default if Git was not found")

    find_package(Git QUIET)
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
            OUTPUT_VARIABLE CEU_CM_GIT_COMMIT_MESSAGE
            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
        if (NOT CEU_CM_GIT_COMMIT_HASH OR NOT CEU_CM_GIT_COMMIT_DATE OR NOT CEU_CM_GIT_COMMIT_MESSAGE)
            message(WARNING "Git found but we're not in a Git repo")
            set(CEU_CM_GIT_COMMIT_HASH
                "UNKNOWN"
                CACHE STRING "Git found but we're not in a Git repo")
            set(CEU_CM_GIT_COMMIT_DATE
                "UNKNOWN"
                CACHE STRING "Git found but we're not in a Git repo")
            set(CEU_CM_GIT_COMMIT_MESSAGE
                "UNKNOWN"
                CACHE STRING "Git found but we're not in a Git repo")
        endif()
    else()
        message(WARNING "Git executable not found.")
    endif()
    message(STATUS "/------------------- Finding and assessing git -------------------\\")
    message(STATUS "|git executable found at ${GIT_EXECUTABLE}")
    message(
        STATUS "|On Git commit ${CEU_CM_GIT_COMMIT_MESSAGE} (${CEU_CM_GIT_COMMIT_HASH}) at ${CEU_CM_GIT_COMMIT_DATE}")
    message(STATUS "\\------------------- Finding and assessing git -------------------/")
endif()
