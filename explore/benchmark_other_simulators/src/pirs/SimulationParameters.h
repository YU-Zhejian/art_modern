#ifndef _SIMULATION_PARAMETERS_H
#define _SIMULATION_PARAMETERS_H

#include <string>
#include "util.h"
#include "InputStream.h"
#include "OutputStream.h"
#include "MaskQvalsByEamss.h"

using std::string;

enum SubstitutionErrorAlgorithm {
	ALGO_DIST,
	ALGO_QTRANS,
};

class SimulationParameters {
public:
	int read_len;
	double coverage;
	double insert_len_mean;
	double insert_len_sd;
	bool jumping;
	bool diploid;
	const char *base_calling_profile_filename;
	const char *gc_bias_profile_filename;
	const char *indel_profile_filename;
	double error_rate;
	enum SubstitutionErrorAlgorithm subst_error_algo;
	casava::demultiplex::QualityMaskMode quality_mask_mode;
	int quality_shift;
	bool simulate_quality_values;
	bool simulate_substitution_errors;
	bool simulate_indels;
	bool simulate_gc_bias;
	bool write_log_files;
	bool user_specified_random_seed;
	uint64_t random_seed;
	int num_simulator_threads;
	string output_directory;
	const char **input_refs;
	unsigned num_input_refs;
	string indiv_name;

	SimulationParameters(int argc, char **argv);
};

class SimulationFiles {
public:
	InputStream  *in_files;
	OutputStream  out_file_1;
	OutputStream  out_file_2;
	OutputStream  insert_distr_log_file;
	OutputStream  info_log_file;
	OutputStream  error_distr_log_file;

	SimulationFiles(const SimulationParameters &params);

	~SimulationFiles() {
		delete[] in_files;
	}
};

class GCBiasProfile;
class BaseCallingProfile;
class IndelProfile;


#endif /* _SIMULATION_PARAMETERS_H */

