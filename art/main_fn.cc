#include "ArtJobPool.hh"
#include <boost/log/trivial.hpp>

#include "ArtConstants.hh"
#include "fasta/FastaStreamBatcher.hh"
#include "fasta/Pbsim3TranscriptBatcher.hh"
#include "main_fn.hh"
#include <fstream>

using namespace std;

namespace labw {
namespace art_modern {

    void print_banner()
    {
        BOOST_LOG_TRIVIAL(info) << "YuZJ Modified ART_Illumina (" << ART_PROGRAM_NAME << ")";
        BOOST_LOG_TRIVIAL(info) << "Based on: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)";
        BOOST_LOG_TRIVIAL(info) << "Originally written by: Weichun Huang <whduke@gmail.com>";
        BOOST_LOG_TRIVIAL(info) << "Modified by: YU Zhejian <Zhejian.23@intl.zju.edu.cn>";
    }

    void generate_all(const ArtParams& art_params)
    {
        ArtJobPool job_pool(art_params);
        if (art_params.art_simulation_mode == SIMULATION_MODE::WGS) {
            // Coverage-based parallelism
            auto coverage_info = art_params.coverage_info.div(art_params.parallel);
            for (int i = 0; i < art_params.parallel; ++i) {
                ArtJobExecutor aje(SimulationJob(art_params.fasta_fetch, coverage_info), art_params);
                job_pool.add(aje);
            }
        } else {
            // Batch-based parallelism
            if (art_params.art_input_file_type == INPUT_FILE_TYPE::FASTA) {
                auto const& coverage_info = art_params.coverage_info;
                if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                    InMemoryFastaStreamBatcher fsb(
                        static_cast<int>(art_params.fasta_fetch->num_seqs() / art_params.parallel + 1),
                        art_params.fasta_fetch);
                    while (true) {
                        auto fa_view = fsb.fetch();

                        if (fa_view.num_seqs() == 0) {
                            break;
                        }
                        ArtJobExecutor aje(
                            SimulationJob(std::make_shared<InMemoryFastaFetch>(fa_view), coverage_info), art_params);
                        job_pool.add(aje);
                    }
                } else {
                    std::ifstream fasta_stream(art_params.input_file_name);
                    FastaStreamBatcher fsb(1000, fasta_stream);
                    while (true) {
                        auto fa_view = fsb.fetch();

                        if (fa_view.num_seqs() == 0) {
                            break;
                        }
                        ArtJobExecutor aje(
                            SimulationJob(std::make_shared<InMemoryFastaFetch>(fa_view), coverage_info), art_params);
                        job_pool.add(aje);
                    }
                    fasta_stream.close();
                }
            }
            if (art_params.art_input_file_type == INPUT_FILE_TYPE::PBSIM3_TEMPLATE) {
                if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                    std::ifstream pbsim3_template_stream(art_params.input_file_name);
                    Pbsim3TranscriptBatcher fsb(1000, pbsim3_template_stream);
                    while (true) {
                        auto fa_view = fsb.fetch();
                        if (fa_view.first.num_seqs() == 0) {
                            break;
                        }
                        ArtJobExecutor aje(
                            SimulationJob(std::make_shared<InMemoryFastaFetch>(fa_view.first), fa_view.second),
                            art_params);
                        job_pool.add(aje);
                    }
                    pbsim3_template_stream.close();
                } else {
                    std::ifstream pbsim3_template_stream(art_params.input_file_name);
                    Pbsim3TranscriptBatcher fsb(1000, pbsim3_template_stream);
                    while (true) {
                        auto fa_view = fsb.fetch();
                        if (fa_view.first.num_seqs() == 0) {
                            break;
                        }
                        ArtJobExecutor aje(
                            SimulationJob(std::make_shared<InMemoryFastaFetch>(fa_view.first), fa_view.second),
                            art_params);
                        job_pool.add(aje);
                    }
                    pbsim3_template_stream.close();
                }
            }
        }
        job_pool.stop();
    }
} // namespace art_modern
} // namespace labw
