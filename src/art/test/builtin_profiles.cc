/**
 * @brief Test whether all the built-in profiles are working.
 *
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

#include "art/builtin_profiles.h"
#include "art/lib/Empdist.hh"

#include "libam_support/Constants.hh"

#include <boost/log/trivial.hpp>

using namespace labw::art_modern;

int main()
{
    for (const bool sep_flag : { true, false }) {
        for (int i = 0; i < N_BUILTIN_PROFILE; i++) {
            const auto is_pe = ENCODED_BUILTIN_PROFILES[i][1][0] != '\0';
            const auto& builtin_profile_name = BUILTIN_PROFILE_NAMES[i];
            BOOST_LOG_TRIVIAL(info) << "Loading builtin profile " << i << ": " << builtin_profile_name;
            auto qdist = Empdist(builtin_profile_name,sep_flag, is_pe);
            qdist.shift_all_emp(0, 0, MIN_QUAL, MAX_QUAL);
            qdist.index();
        }
    }
}
