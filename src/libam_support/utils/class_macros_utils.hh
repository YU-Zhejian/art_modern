#pragma once
// NOLINTBEGIN

#define DELETE_COPY(classname)                                                                                         \
    classname(const classname&) = delete;                                                                              \
    classname& operator=(const classname&) = delete;
#define DELETE_MOVE(classname)                                                                                         \
    classname(classname&&) = delete;                                                                                   \
    classname& operator=(classname&&) = delete;

#define DEFAULT_COPY(classname)                                                                                        \
    classname(const classname&) = default;                                                                             \
    classname& operator=(const classname&) = default;
#define DEFAULT_MOVE(classname)                                                                                        \
    classname(classname&&) = default;                                                                                  \
    classname& operator=(classname&&) = default;
// NOLINTEND
