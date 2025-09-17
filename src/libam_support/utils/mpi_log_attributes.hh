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
#include <boost/log/attributes/attribute.hpp>
#include <boost/log/attributes/attribute_cast.hpp>

namespace labw::art_modern {

/**
 * An integer rank, returning MPI_UNAVAILABLE_RANK if unavailable.
 */
class MPIRankLoggerAttribute : public boost::log::attribute {
public:
    MPIRankLoggerAttribute();
    [[maybe_unused]] explicit MPIRankLoggerAttribute(boost::log::attributes::cast_source const& source);
};

/**
 * A string representing the hostname of the MPI process. Reyturn "N/A" if unavailable.
 */
class MPIHostNameLoggerAttribute : public boost::log::attribute {
public:
    MPIHostNameLoggerAttribute();
    [[maybe_unused]] explicit MPIHostNameLoggerAttribute(boost::log::attributes::cast_source const& source);
};
} // namespace labw::art_modern
