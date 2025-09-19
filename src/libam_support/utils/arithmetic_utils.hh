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

namespace labw::art_modern {
/** @brief Returns the maximum of two values.
 *
 * This function takes two values of the same type and returns the larger of the two.
 *
 * @tparam T The type of the input values. Must support comparison operators.
 * @param a The first value to compare.
 * @param b The second value to compare.
 * @return The maximum of the two input values.
 */
template <typename T> T am_max(T a, T b)
{
    if (a > b) {
        return a;
    }
    return b;
}
/** @brief Returns the minimum of two values.
 *
 * This function takes two values of the same type and returns the smaller of the two.
 *
 * @tparam T The type of the input values. Must support comparison operators.
 * @param a The first value to compare.
 * @param b The second value to compare.
 * @return The minimum of the two input values.
 */
template <typename T> T am_min(T a, T b)
{
    if (a < b) {
        return a;
    }
    return b;
}

} // namespace labw::art_modern
