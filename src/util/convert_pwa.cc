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

#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>

#include "libam_support/Constants.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/OutputDispatcher.hh"
#include "libam_support/ref/batcher/Pbsim3TranscriptBatcher.hh"
#include "libam_support/ref/fetch/FaidxFetch.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/utils/fs_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/seq_utils.hh"

namespace po = boost::program_options;
using namespace labw::art_modern;

const char ARG_INPUT_PWA_FILE_NAME[] = "i-pwa";
const char ARG_INPUT_REF_FILE_NAME[] = "i-ref";
const char ARG_INPUT_FILE_PARSER[] = "i-parser";
const char ARG_INPUT_FILE_TYPE[] = "i-type";
const char ARG_FCOV[] = "i-fcov";
const char ARG_BATCH_SIZE[] = "i-batch_size";
const char ARG_IS_PE[] = "paired-end";

INPUT_FILE_TYPE get_input_file_type(const std::string& input_file_type_str, const std::string& input_file_name)
{
    if (input_file_type_str == INPUT_FILE_TYPE_FASTA) {
        return INPUT_FILE_TYPE::FASTA;
    } else if (input_file_type_str == INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS) {
        return INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS;
    } else if (input_file_type_str == INPUT_FILE_TYPE_AUTO) {
        for (const auto& fasta_file_end : std::vector<std::string> { "fna", "fsa", "fa", "fasta" }) {
            if (ends_with(input_file_name, fasta_file_end)) {
                return INPUT_FILE_TYPE::FASTA;
            }
        }
        BOOST_LOG_TRIVIAL(fatal) << "Automatic inference of input file type failed! Modify value of this param (--"
                                 << ARG_INPUT_FILE_TYPE << ") to be one of " << INPUT_FILE_TYPE_FASTA << ", "
                                 << INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS << ".";
        abort_mpi();
    } else {
        BOOST_LOG_TRIVIAL(fatal) << "Input file type (--" << ARG_INPUT_FILE_TYPE << ") should be one of "
                                 << INPUT_FILE_TYPE_FASTA << ", " << INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS << ", "
                                 << INPUT_FILE_TYPE_AUTO << ".";
        abort_mpi();
    }
}

INPUT_FILE_PARSER get_input_file_parser(const std::string& input_file_parser_str)
{
    if (input_file_parser_str == INPUT_FILE_PARSER_MEMORY) {
        return INPUT_FILE_PARSER::MEMORY;
    } else if (input_file_parser_str == INPUT_FILE_PARSER_HTSLIB) {
        return INPUT_FILE_PARSER::HTSLIB;
    } else if (input_file_parser_str == INPUT_FILE_PARSER_STREAM) {
        return INPUT_FILE_PARSER::STREAM;
    } else {
        BOOST_LOG_TRIVIAL(fatal) << "Input file parser (--" << ARG_INPUT_FILE_PARSER << ") should be one of "
                                 << INPUT_FILE_PARSER_MEMORY << ", " << INPUT_FILE_PARSER_HTSLIB << ", "
                                 << INPUT_FILE_PARSER_STREAM << ".";
        abort_mpi();
    }
}

po::options_description get_po_desc(const OutputDispatcherFactory& factory)
{
    po::options_description required_opts("Required Options");

    required_opts.add_options()(ARG_INPUT_PWA_FILE_NAME, po::value<std::string>(), "input PWA file");
    required_opts.add_options()(
        ARG_IS_PE, "whether the input PWA was generated using paired-end (PE) or mate-pair (MP) sequencing");
    required_opts.add_options()(ARG_INPUT_FILE_PARSER, po::value<std::string>()->default_value(INPUT_FILE_PARSER_AUTO),
        (std::string() + "input file parser, should be " + INPUT_FILE_PARSER_MEMORY + ", " + INPUT_FILE_PARSER_HTSLIB
            + ", " + INPUT_FILE_PARSER_STREAM + ".")
            .c_str());
    required_opts.add_options()(ARG_INPUT_FILE_TYPE, po::value<std::string>()->default_value(INPUT_FILE_TYPE_AUTO),
        (std::string() + "input file type, should be " + INPUT_FILE_TYPE_AUTO + ", " + INPUT_FILE_TYPE_FASTA + ", "
            + INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS + ".")
            .c_str());

    required_opts.add_options()(ARG_INPUT_REF_FILE_NAME, po::value<std::string>(),
        "the filename of input reference genome, reference "
        "transcriptome, or templates");

    auto po_desc = po::options_description("Options");
    po_desc.add(required_opts);
    factory.patch_options(po_desc);
    return po_desc;
}

int main(int argc, char** argv)
{
    const std::vector<std::string> args { argv, argv + argc };
    po::variables_map vm;
    po::options_description required_opts("Required Options");

    const auto& factory = get_output_dispatcher_factory();
    const auto& po_desc = get_po_desc(factory);

    try {
        store(po::parse_command_line(argc, argv, po_desc), vm);
        notify(vm);
    } catch (const std::exception& exp) {
        BOOST_LOG_TRIVIAL(fatal) << exp.what();
        std::cout << po_desc << std::endl;
        abort_mpi();
    }

    BaseReadOutput* out;

    const auto& input_pwa_file = vm[ARG_INPUT_PWA_FILE_NAME].as<std::string>();
    validate_input_filename(input_pwa_file, ARG_INPUT_PWA_FILE_NAME);

    const auto& input_ref_file = vm[ARG_INPUT_REF_FILE_NAME].as<std::string>();
    validate_input_filename(input_ref_file, ARG_INPUT_REF_FILE_NAME);

    const auto input_file_type = get_input_file_type(vm[ARG_INPUT_FILE_TYPE].as<std::string>(), input_ref_file);
    const auto input_file_parser = get_input_file_parser(vm[ARG_INPUT_FILE_PARSER].as<std::string>());
    if (input_file_parser == INPUT_FILE_PARSER::STREAM) {
        InMemoryFastaFetch fetch;
        out = factory.create(vm, &fetch, args);
    } else if (input_file_type == INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS) {
        if (input_file_parser == INPUT_FILE_PARSER::MEMORY) {
            std::ifstream input_stream(input_ref_file);
            Pbsim3TranscriptBatcher pbsim3_transcript_batcher(std::numeric_limits<int>::max(), input_stream);
            const auto& [fetch_, cov] = pbsim3_transcript_batcher.fetch();
            out = factory.create(vm, &fetch_, args);
            input_stream.close();
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "Input file type PBSIM3_TRANSCRIPTS cannot be parsed using HTSLib.";
            abort_mpi();
        }
    } else {
        BaseFastaFetch* fetch;
        if (input_file_parser == INPUT_FILE_PARSER::MEMORY) {
            fetch = new InMemoryFastaFetch(input_ref_file);
        } else {
            fetch = new FaidxFetch(input_ref_file);
        }
        out = factory.create(vm, fetch, args);
        delete fetch;
    }

    const bool is_pe = vm.count(ARG_IS_PE) > 0;

    std::ifstream ifs(input_pwa_file);
    std::string line;
    std::array<std::string, PairwiseAlignment::NUM_LINES> serialized;
    while (!ifs.eof()) {
        std::getline(ifs, line);
        if (line.empty()) {
            break;
        }
        serialized[0] = std::move(line);
        for (auto i = 1; i < PairwiseAlignment::NUM_LINES; ++i) {
            std::getline(ifs, line);
            serialized[i] = std::move(line);
        }
        if (!is_pe) {
            out->writeSE(PairwiseAlignment::deserialize(serialized));
        } else {
            const PairwiseAlignment& pwa1 = PairwiseAlignment::deserialize(serialized);
            for (auto i = 0; i < PairwiseAlignment::NUM_LINES; ++i) {
                std::getline(ifs, line);
                serialized[i] = std::move(line);
            }
            const PairwiseAlignment& pwa2 = PairwiseAlignment::deserialize(serialized);
            out->writePE(pwa1, pwa2);
        }
    }
    ifs.close();
    out->close();
    delete out;
    exit_mpi(EXIT_SUCCESS);
}
