#pragma once
// NOLINTBEGIN

#define DELETE_COPY(classname)                                                                                         \
    classname(const classname&) = delete;                                                                              \
    classname& operator=(const classname&) = delete;
#define DELETE_MOVE(classname)                                                                                         \
    classname(classname&&) = delete;                                                                                   \
    classname& operator=(classname&&) = delete;

// NOLINTEND
