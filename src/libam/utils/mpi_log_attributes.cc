#include "libam/utils/mpi_log_attributes.hh"

#include "libam/Constants.hh"

#include <boost/log/attributes/attribute.hpp>
#include <boost/log/attributes/attribute_cast.hpp>
#include <boost/log/attributes/attribute_value.hpp>
#include <boost/log/attributes/attribute_value_impl.hpp>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <string>

namespace labw::art_modern {

class MPIRankLoggerAttributeImpl final : public boost::log::attribute::impl {
public:
    boost::log::attribute_value get_value() override
    {
#ifdef WITH_MPI
        int mpi_finalized_flag;
        MPI_Finalized(&mpi_finalized_flag);
        if (mpi_finalized_flag) {
            return boost::log::attributes::make_attribute_value(MPI_UNAVAILABLE_RANK);
        }
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return boost::log::attributes::make_attribute_value(rank);
#else
        return boost::log::attributes::make_attribute_value(MPI_UNAVAILABLE_RANK);
#endif
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
#ifdef WITH_MPI
        int mpi_finalized_flag;
        MPI_Finalized(&mpi_finalized_flag);
        if (mpi_finalized_flag) {
            return boost::log::attributes::make_attribute_value(std::string("N/A"));
        }
        char hostname[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(hostname, &name_len);
        return boost::log::attributes::make_attribute_value(std::string(hostname, name_len));
#else
        return boost::log::attributes::make_attribute_value(std::string("N/A"));
#endif
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