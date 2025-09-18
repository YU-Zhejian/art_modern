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

#include "art_modern_config.h" // NOLINT: For WITH_PROTOBUF

#include "libam_support/utils/protobuf_utils.hh"

#ifdef WITH_PROTOBUF
#include <google/protobuf/stubs/common.h>
#endif

void labw::art_modern::validate_protobuf_version()
{
#ifdef WITH_PROTOBUF
    GOOGLE_PROTOBUF_VERIFY_VERSION;
#endif
}
