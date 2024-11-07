#include "utils/mpi_log_attributes.hh"

#include <boost/log/attributes/attribute_value_impl.hpp>
#include <mpi.h>

namespace labw::art_modern {

class MPIRankLoggerAttributeImpl : public boost::log::attribute::impl {
public:
    boost::log::attribute_value get_value() override
    {
        int mpi_finalized_flag;
        MPI_Finalized(&mpi_finalized_flag);
        if (mpi_finalized_flag) {
            return boost::log::attributes::make_attribute_value(MPI_UNAVAILABLE_RANK);
        }
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return boost::log::attributes::make_attribute_value(rank);
    }
};

[[maybe_unused]] MPIRankLoggerAttribute::MPIRankLoggerAttribute(const boost::log::attributes::cast_source& source)
    : boost::log::attribute(source.as<MPIRankLoggerAttributeImpl>())
{
}
MPIRankLoggerAttribute::MPIRankLoggerAttribute()
    : boost::log::attribute(new MPIRankLoggerAttributeImpl())
{
}

class MPIHostNameLoggerAttributeImpl : public boost::log::attribute::impl {
public:
    boost::log::attribute_value get_value() override
    {
        int mpi_finalized_flag;
        MPI_Finalized(&mpi_finalized_flag);
        if (mpi_finalized_flag) {
            return boost::log::attributes::make_attribute_value(std::string("N/A"));
        }
        char hostname[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(hostname, &name_len);
        return boost::log::attributes::make_attribute_value(std::string(hostname, name_len));
    }
};

MPIHostNameLoggerAttribute::MPIHostNameLoggerAttribute(const boost::log::attributes::cast_source& source)
    : boost::log::attribute(source.as<MPIHostNameLoggerAttributeImpl>())
{
}
MPIHostNameLoggerAttribute::MPIHostNameLoggerAttribute()
    : boost::log::attribute(new MPIHostNameLoggerAttributeImpl())
{
}
}