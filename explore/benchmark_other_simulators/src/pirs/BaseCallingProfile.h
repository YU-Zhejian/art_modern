#ifndef _BASE_CALLING_PROFILE_H
#define _BASE_CALLING_PROFILE_H

#include <inttypes.h>

#include "Profile.h"
#include "SimulationParameters.h"

class Random;

/* Profile for simulating substitution errors and quality values in the reads.
 * */
class BaseCallingProfile : public Profile {
private:
	void load_dimensions_header();
	friend void load_dist_matrix_line(char *line, Profile &_profile);
	friend void load_qtrans_matrix_line(char *line, Profile &_profile);
	int get_cycle(const char *cycle_str);
public:
	int num_ref_bases;
	int num_cycles;
	int num_seq_bases;
	int num_quality_scores;
	int read_len;
	double statistical_average_error_rate;
	uint64_t ***dist_matrix;
	uint64_t ****dist_called_base_matrix;
	uint64_t ***qtrans_matrix;

	BaseCallingProfile(const SimulationParameters &params);
	~BaseCallingProfile();
	char call(int cycle, char raw_base, int prev_qscore,
		  int &qscore, Random &rgen) const;

};

#define DEFAULT_BASE_CALLING_PROFILE PKGDATADIR "/Base-Calling_Profiles/humNew.PE100.matrix.gz"

#endif /* ifndef _BASE_CALLING_PROFILE_H */
