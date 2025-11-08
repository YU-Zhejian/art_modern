/**
 * Copyright 2008-2016 Weichun Huang <whduke@gmail.com>
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

#include "art/lib/ArtJobExecutor.hh"

#include "art/lib/ArtConstants.hh"
#include "art/lib/ArtContig.hh"
#include "art/lib/ArtRead.hh"
#include "art/lib/Rprob.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/out/OutputDispatcher.hh"
#include "libam_support/utils/arithmetic_utils.hh"
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/depth_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/si_utils.hh"

#include <fmt/format.h>

#include <boost/log/trivial.hpp>

#include <htslib/hts.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <thread>
#include <utility>

namespace labw::art_modern {
namespace {

    class AJEReporter {
    public:
        explicit AJEReporter(ArtJobExecutor& aje, const std::size_t reporting_interval_seconds)
            : aje_(aje)
            , reporting_interval_seconds_(reporting_interval_seconds)
        {
        }

        DELETE_COPY(AJEReporter)
        DELETE_MOVE(AJEReporter)
        ~AJEReporter() { stop(); }

        void stop()
        {
            should_stop_ = true;
            if (thread_.joinable()) {
                thread_.join();
            }
        }
        void start() { thread_ = std::thread(&AJEReporter::job_, this); }

    private:
        void job_() const
        {
            std::this_thread::sleep_for(std::chrono::seconds(reporting_interval_seconds_));
            while (!should_stop_) {
                BOOST_LOG_TRIVIAL(info) << "AJEReporter: Job " << aje_.thread_info();
                std::this_thread::sleep_for(std::chrono::seconds(reporting_interval_seconds_));
            }
        }

        ArtJobExecutor& aje_;
        std::atomic<bool> should_stop_ { false };
        std::thread thread_;
        std::size_t reporting_interval_seconds_;
    };

} // namespace
void ArtJobExecutor::generate_(const am_readnum_t targeted_num_reads, const bool is_positive, ArtContig& art_contig,
    am_readnum_t& read_id, Rprob& rprob_)
{
    current_contig_name_ = art_contig.seq_name;
    current_n_fails_ = 0;
    current_n_reads_generated_ = 0;
    current_max_tolerence_ = am_max(static_cast<am_readnum_t>(5),
        static_cast<am_readnum_t>(static_cast<double>(targeted_num_reads) * MAX_TRIAL_RATIO_BEFORE_FAIL));
    bool retv = false;
    current_n_reads_left_ = targeted_num_reads;
    while (current_n_reads_left_ > 0) {
        if (art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::SE) {
            retv = generate_se_(art_contig, is_positive, read_id, rprob_);
        } else {
            retv = generate_pe_(art_contig, is_positive, read_id, rprob_);
        }
        if (retv) {
            current_n_reads_left_ -= num_reads_to_reduce_;
            total_num_reads_generated_ += num_reads_to_reduce_;
            current_n_reads_generated_ += num_reads_to_reduce_;
            read_id++;
        } else {
            current_n_fails_++;
        }
        if (current_n_fails_ >= current_max_tolerence_) {
            BOOST_LOG_TRIVIAL(debug) << "Failed to generate reads for " << art_contig.seq_name << " sized "
                                     << art_contig.seq_size << " after " << current_max_tolerence_ << " attempts.";
            break;
        }
    }
}

bool ArtJobExecutor::generate_pe_(
    ArtContig& art_contig, const bool is_plus_strand, const am_readnum_t current_num_reads, Rprob& rprob_)
{
    const std::string read_name = fmt::format(
        "{}:{}:{}:{}:{}", art_contig.seq_name, art_params_.id, job_.job_id, mpi_rank_str_, current_num_reads);

    ArtRead read_1(art_params_, art_contig.seq_name, read_name, true, rprob_);
    ArtRead read_2(art_params_, art_contig.seq_name, read_name, false, rprob_);
    art_contig.generate_read_pe(is_plus_strand, read_1, read_2);
    read_1.generate_snv_on_qual();
    read_2.generate_snv_on_qual();
    if (require_alignment_) {
        read_1.generate_pairwise_aln();
        read_2.generate_pairwise_aln();
    }
    if (!(read_1.is_good() && read_2.is_good())) {
        return false;
    }

    output_dispatcher_->writePE(token_ring_, read_1.to_pwa(), read_2.to_pwa());

    return true;
}

bool ArtJobExecutor::generate_se_(
    ArtContig& art_contig, const bool is_plus_strand, const am_readnum_t current_num_reads, Rprob& rprob_)
{
    auto read_id = fmt::format(
        "{}:{}:{}:{}:{}", art_contig.seq_name, art_params_.id, job_.job_id, mpi_rank_str_, current_num_reads);
    ArtRead art_read(art_params_, art_contig.seq_name, std::move(read_id), true, rprob_);
    art_contig.generate_read_se(is_plus_strand, art_read);
    art_read.generate_snv_on_qual();
    if (require_alignment_) {
        art_read.generate_pairwise_aln();
    }
    if (!art_read.is_good()) {
        return false;
    }
    output_dispatcher_->writeSE(token_ring_, art_read.to_pwa());
    return true;
}

ArtJobExecutor::ArtJobExecutor(
    SimulationJob&& job, const ArtParams& art_params, const std::shared_ptr<OutputDispatcher>& output_dispatcher)
    : art_params_(art_params)
    , job_(std::move(job))
    , mpi_rank_str_(mpi_rank_s())
    , output_dispatcher_(output_dispatcher)
    , num_reads_to_reduce_(art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::SE ? 1 : 2)
    , require_alignment_(output_dispatcher->require_alignment())
    , token_ring_(output_dispatcher->get_producer_tokens())
{
}

bool ArtJobExecutor::is_running() const { return is_running_; }

void ArtJobExecutor::operator()()
{
    is_running_ = true;
    const auto num_contigs = job_.fasta_fetch->num_seqs();
    if (num_contigs == 0) {
        is_running_ = false;
        return;
    }

    AJEReporter reporter(*this, art_params_.art_job_executor_reporting_interval_seconds);
    reporter.start();
    hts_pos_t accumulated_contig_len = 0;

    Rprob rprob_(art_params_.pe_frag_dist_mean, art_params_.pe_frag_dist_std_dev, art_params_.read_len_1,
        art_params_.read_len_2);

    BOOST_LOG_TRIVIAL(info) << "Starting simulation for job " << job_.job_id << " with " << num_contigs << " contigs";
    for (decltype(job_.fasta_fetch->num_seqs()) seq_id = 0; seq_id < num_contigs; ++seq_id) {
        const auto& contig_name = job_.fasta_fetch->seq_name(seq_id);
        const auto contig_size = job_.fasta_fetch->seq_len(seq_id);
        if (contig_size < art_params_.read_len_max || /** unlikely */ contig_size == 0) {
            BOOST_LOG_TRIVIAL(debug) << "the reference sequence " << contig_name << " (length "
                                     << job_.fasta_fetch->seq_len(seq_id)
                                     << "bps ) is skipped as it < the defined maximum read length ("
                                     << art_params_.read_len_max << " bps)";
            continue;
        }
        accumulated_contig_len += contig_size;

        ArtContig art_contig(job_.fasta_fetch, seq_id, art_params_, rprob_);

        am_readnum_t num_pos_reads = 0;
        am_readnum_t num_neg_reads = 0;
        const auto cov_pos = job_.coverage_info->coverage_positive(contig_name);
        const auto cov_neg = job_.coverage_info->coverage_negative(contig_name);
        if (art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE) {
            const auto [npr, nnr] = calculate_num_reads_template(cov_pos, cov_neg, num_reads_to_reduce_);
            num_pos_reads = npr;
            num_neg_reads = nnr;
        } else {
            const auto [npr, nnr]
                = calculate_num_reads_old(contig_size, art_params_.read_len_1, cov_pos, cov_neg, num_reads_to_reduce_);
            num_pos_reads = npr;
            num_neg_reads = nnr;
        }

        if (num_pos_reads + num_neg_reads == 0) {
            BOOST_LOG_TRIVIAL(debug) << "No read will be generated for the reference sequence " << contig_name
                                     << " due to insufficient coverage (pos=" << cov_pos << ", neg=" << cov_neg << ")";
            continue;
        }
        am_readnum_t read_id = 0;
        generate_(num_pos_reads, true, art_contig, read_id, rprob_);
        generate_(num_neg_reads, false, art_contig, read_id, rprob_);
    }
    if (accumulated_contig_len == 0) {
        // Happens when all contigs are shorter than read length
        BOOST_LOG_TRIVIAL(info) << "Finished simulation for job " << job_.job_id
                                << " with N/A reads (mean depth=N/A) generated.";
    } else {
        // Calculate real coverage
        double real_coverage = 0;
        if (art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE) {
            // TODO: Still needs checking
            real_coverage = total_num_reads_generated_;
        } else {
            if (art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::SE) {
                real_coverage = static_cast<double>(total_num_reads_generated_) * art_params_.read_len_1
                    / static_cast<double>(accumulated_contig_len);
            } else {
                real_coverage = static_cast<double>(total_num_reads_generated_)
                    * (art_params_.read_len_1 + art_params_.read_len_2) / static_cast<double>(accumulated_contig_len)
                    / 2.0;
            }
        }
        BOOST_LOG_TRIVIAL(info) << "Finished simulation for job " << job_.job_id << " with "
                                << to_si(total_num_reads_generated_) << " reads (mean depth=" << real_coverage
                                << ") generated.";
    }
    reporter.stop();
    is_running_ = false;
}
ArtJobExecutor::ArtJobExecutor(ArtJobExecutor&& other) noexcept
    : art_params_(other.art_params_)
    , job_(std::move(other.job_))
    , mpi_rank_str_(other.mpi_rank_str_)
    , output_dispatcher_(std::move(other.output_dispatcher_))
    , num_reads_to_reduce_(other.num_reads_to_reduce_)
    , require_alignment_(other.require_alignment_)
    , token_ring_(std::move(other.token_ring_))
{
}
std::string ArtJobExecutor::thread_info() const
{
    return fmt::format("{}:{} | ON: '{}' | SUCCESS: current={}, left={} | FAIL: current={}, max={} | TOTAL: {}",
        job_.job_id, mpi_rank_str_, current_contig_name_, to_si(current_n_reads_generated_),
        to_si(current_n_reads_left_), to_si(current_n_fails_), to_si(current_max_tolerence_),
        to_si(total_num_reads_generated_));
}

} // namespace labw::art_modern
