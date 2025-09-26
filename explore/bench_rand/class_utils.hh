#pragma once

// Define a macro that deletes copy and move constructors and assignment operators
#define DELETE_COPY_MOVE(ClassName)                                                                                    \
    ClassName(const ClassName&) = delete;                                                                              \
    ClassName(ClassName&&) = delete;                                                                                   \
    ClassName& operator=(const ClassName&) = delete;                                                                   \
    ClassName& operator=(ClassName&&) = delete;
