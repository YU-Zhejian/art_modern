/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#pragma once
// NOLINTBEGIN

#define DELETE_COPY_ASSIGNMENT(classname) classname& operator=(const classname&) = delete;
#define DELETE_MOVE_ASSIGNMENT(classname) classname& operator=(classname&&) = delete;

#define DELETE_COPY_CONSTRUCTOR(classname) classname(const classname&) = delete;
#define DELETE_MOVE_CONSTRUCTOR(classname) classname(classname&&) = delete;

/**
 * @brief Macro to delete copy constructor and copy assignment operator for a class.
 * @param classname  The name of the class for which to delete copy operations.
 */
#define DELETE_COPY(classname)                                                                                         \
    DELETE_COPY_ASSIGNMENT(classname)                                                                                  \
    DELETE_COPY_CONSTRUCTOR(classname)
/**
 * @brief  Macro to delete move constructor and move assignment operator for a class.
 * @param classname  The name of the class for which to delete move operations.
 */
#define DELETE_MOVE(classname)                                                                                         \
    DELETE_MOVE_ASSIGNMENT(classname)                                                                                  \
    DELETE_MOVE_CONSTRUCTOR(classname)

/**
 *  @brief Macro to default copy constructor and copy assignment operator for a class.
 * @param classname  The name of the class for which to default copy operations.
 */
#define DEFAULT_COPY(classname)                                                                                        \
    classname(const classname&) = default;                                                                             \
    classname& operator=(const classname&) = default;
/**
 *  @brief Macro to default move constructor and move assignment operator for a class.
 * @param classname  The name of the class for which to default move operations.
 */
#define DEFAULT_MOVE(classname)                                                                                        \
    classname(classname&&) = default;                                                                                  \
    classname& operator=(classname&&) = default;

#define DEFAULT_DESTRUCTOR(classname) ~classname() = default;
#define DELETE_DESTRUCTOR(classname) ~classname() = delete;
// NOLINTEND
