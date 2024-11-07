#pragma once

#include <boost/log/attributes/attribute.hpp>
#include <boost/log/attributes/attribute_cast.hpp>
namespace labw::art_modern {

class MPIRankLoggerAttribute : public boost::log::attribute {
public:
    MPIRankLoggerAttribute();
    explicit MPIRankLoggerAttribute(boost::log::attributes::cast_source const& source);
};
}
