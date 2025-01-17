#include "art_modern_config.h" // NOLINT: For WITH_PROTOBUF

#include "libam/utils/protobuf_utils.hh"

#ifdef WITH_PROTOBUF
#include <google/protobuf/stubs/common.h>
#endif

void labw::art_modern::validate_protobuf_version()
{
#ifdef WITH_PROTOBUF
    GOOGLE_PROTOBUF_VERIFY_VERSION;
#endif
}
