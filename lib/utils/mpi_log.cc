#include "utils/mpi_log.hh"

#include <boost/log/attributes/attribute_value_impl.hpp>
#include <mpi.h>

namespace labw::art_modern {

class MPIRankLoggerAttributeImpl : public boost::log::attribute::impl {
public:
    boost::log::attribute_value get_value() override
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return boost::log::attributes::make_attribute_value(rank);
    }
};

MPIRankLoggerAttribute::MPIRankLoggerAttribute(const boost::log::attributes::cast_source& source)
    : boost::log::attribute(source.as<MPIRankLoggerAttributeImpl>())
{
}
MPIRankLoggerAttribute::MPIRankLoggerAttribute()
    : boost::log::attribute(new MPIRankLoggerAttributeImpl())
{
}
}