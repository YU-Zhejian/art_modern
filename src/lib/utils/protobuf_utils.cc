#include "protobuf_utils.hh"
#include "art_modern_config.h" // For WITH_PROTOBUF

#ifdef WITH_PROTOBUF
#include <google/protobuf/stubs/common.h>
#endif

void labw::art_modern::validate_protobuf_version()
{
#ifdef WITH_PROTOBUF
    GOOGLE_PROTOBUF_VERIFY_VERSION;
#endif
}
