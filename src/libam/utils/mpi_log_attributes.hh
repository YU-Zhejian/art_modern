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
