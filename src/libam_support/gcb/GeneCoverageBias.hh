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

#pragma once

#include <cstdint>


    namespace labw::art_modern {
        enum class GENE_AMPLIFY_MODE: std::uint8_t {
            /** No gene coverage bias */
            NONE,
            /** 10X Genomics 3' scRNA-Seq protocol */
            TEN_TIMES_3P,
            /** 10X Genomics 5' scRNA-Seq protocol */
            TEN_TIMES_5P
        };

        class GeneCoverageBias {
        public:
            /**
             * Default constructor.
             */
            GeneCoverageBias() = default;
        private:
            /**
             * Parameters for 2D-KDE model for gene truncation
             */
            void* gene_trunc_kde_model_ = nullptr;
            GENE_AMPLIFY_MODE gene_amplify_mode_ = GENE_AMPLIFY_MODE::NONE;
        };

    } // namespace labw::art_modern


