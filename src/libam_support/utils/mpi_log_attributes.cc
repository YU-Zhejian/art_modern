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

#include "art_modern_config.h" // NOLINT: For WITH_MPI

#include "libam_support/utils/mpi_log_attributes.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/attributes/attribute.hpp>
#include <boost/log/attributes/attribute_cast.hpp>
#include <boost/log/attributes/attribute_value.hpp>
#include <boost/log/attributes/attribute_value_impl.hpp>

namespace labw::art_modern {

class MPIRankLoggerAttributeImpl final : public boost::log::attribute::impl {
public:
    boost::log::attribute_value get_value() override
    {
        return boost::log::attributes::make_attribute_value(mpi_rank_s());
    }
};

[[maybe_unused]] MPIRankLoggerAttribute::MPIRankLoggerAttribute(const boost::log::attributes::cast_source& source)
    : attribute(source.as<MPIRankLoggerAttributeImpl>())
{
}
MPIRankLoggerAttribute::MPIRankLoggerAttribute()
    : attribute(new MPIRankLoggerAttributeImpl())
{
}

class MPIHostNameLoggerAttributeImpl final : public boost::log::attribute::impl {
public:
    boost::log::attribute_value get_value() override
    {
        return boost::log::attributes::make_attribute_value(mpi_hostname());
    }
};

MPIHostNameLoggerAttribute::MPIHostNameLoggerAttribute(const boost::log::attributes::cast_source& source)
    : attribute(source.as<MPIHostNameLoggerAttributeImpl>())
{
}
MPIHostNameLoggerAttribute::MPIHostNameLoggerAttribute()
    : attribute(new MPIHostNameLoggerAttributeImpl())
{
}
} // namespace labw::art_modern
