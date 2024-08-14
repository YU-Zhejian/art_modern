#pragma once
#include "Empdist.hh"
#include <boost/program_options.hpp>
#include <map>
#include <string>
#include <vector>

namespace labw {
namespace art_modern {

    class ArtParams {
    public:
        ART_SIMULATION_MODE art_simulation_mode;
        ART_LIB_CONST_MODE art_lib_const_mode;
        int min_qual;
        int max_qual;
        bool parallel_on_read;
        int parallel;
        std::vector<std::string> _args;

        /**
         * @brief a negative value means no limit
         */
        int max_indel;
        bool sep_flag = false;
        bool cigar_use_m = false;
        bool no_sam;
        bool stream;
        std::string id;
        std::string qual_file_1;
        std::string qual_file_2;
        std::string seq_file;
        std::string out_file_prefix;
        std::map<std::string, double> sequencing_depth;
        double uniform_sequencing_depth = 0.0;
        std::string p_cigar;
        std::string fcov;
        int read_len;
        int q_shift_1;
        int q_shift_2;
        int pe_frag_dist_mean;
        double pe_frag_dist_std_dev;
        double ins_rate_1 = DEFAULT_INS_RATE_1;
        double del_rate_1 = DEFAULT_DEL_RATE_1;
        double ins_rate_2 = DEFAULT_INS_RATE_2;
        double del_rate_2 = DEFAULT_DEL_RATE_2;
        std::vector<double> per_base_ins_rate_1;
        std::vector<double> per_base_del_rate_1;
        std::vector<double> per_base_ins_rate_2;
        std::vector<double> per_base_del_rate_2;
        double err_prob[HIGHEST_QUAL];
        boost::program_options::options_description po_desc;

        ArtParams();

        void parse_args(const std::vector<std::string>& args);

        void validate_args();

        void shift_emp(std::map<int, int> map_to_process, int q_shift) const;

        Empdist read_emp() const;

        /**
         * Get file name for the FASTQ file or FASTQ file 1 for PE simulation.
         * @param gene_name As described.
         * @return As described.
         */
        std::string fqfile1(const std::string& gene_name) const;

        /**
         * Get file name for the FASTQ file.
         * @param gene_name As described.
         * @return As described.
         */
        std::string fqfile(const std::string& gene_name) const;

        /**
         * Get file name for the FASTQ file 2.
         * @param gene_name As described.
         * @return As described.
         */
        std::string fqfile2(const std::string& gene_name) const;

        /**
         * Get file name for the SAM file.
         * @param gene_name As described.
         * @return As described.
         */
        std::string samfile(const std::string& gene_name) const;

        /**
         * Print parameter parsing results.
         */
        void print_params() const;
        void print_help() const;
    };

} // namespace art_modern
} // namespace labw
