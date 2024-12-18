#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "BaseCallingProfile.h"
#include "Random.h"

static bool get_int(const char *tok, int &n)
{
	char *tmp;
	if (!tok)
		return false;
	n = strtol(tok, &tmp, 10);
	if (tmp == tok)
		return false;
	return true;
}

int BaseCallingProfile::get_cycle(const char *cycle_str)
{
	int cycle;

	if (!get_int(cycle_str, cycle))
		return -1;

	if (cycle >= 1 && cycle <= params.read_len)
		return cycle;

	cycle -= num_cycles / 2;
	if (cycle >= 1 && cycle <= params.read_len)
		return cycle + params.read_len;

	return -1;
}

//
// Load @len 64-bit unsigned integers from the strtok() string into the array
// @row.  Return true on success; false if there are not enough tokens, or if
// there is an invalid token.
//
static bool load_matrix_row(uint64_t row[], size_t len)
{
	char *tmp, *tok;
	while (len--) {
		tok = strtok(NULL, " \t");
		if (!tok)
			return false;
		*row++ = strtoull(tok, &tmp, 10);
		if (tmp == tok)
			return false;
	}
	return true;
}

//
// Converts quality score to the corresponding error rate.
//
static inline double qscore_to_error_rate(int qscore)
{
	return pow(10.0, -qscore / 10.0);
}

//
// Converts error rate to the corresponding quality score.
//
static inline int error_rate_to_qscore(double error_rate)
{
	return int(-10 * log10(error_rate) + 0.5);
}

//
// Builds an array @qval_to_qval that maps quality scores to quality scores,
// such that the mapped-to quality scores indicate an error rate @factor times
// higher than the original quality score.
//
static void qscore_calibrate(double factor, unsigned qval_to_qval[],
			     unsigned num_quality_scores)
{
	for (unsigned i = 0; i < num_quality_scores; i++) {

		int new_qscore = error_rate_to_qscore(factor * qscore_to_error_rate(i));
		
		if(new_qscore < 2)
			qval_to_qval[i] = 2;
		else if (new_qscore >= (int)num_quality_scores)
			qval_to_qval[i] = num_quality_scores - 1;
		else
			qval_to_qval[i] = new_qscore;
	}
}

//
// Builds an array @qval_to_qval that re-maps quality scores so that the average
// error rate of the quality score distribution @qscore_ratio_distr becomes
// close to @desired_error_rate.
//
static void compute_qval_to_qval(unsigned qval_to_qval[],
				 const double qscore_ratio_distr[],
				 unsigned num_quality_scores,
				 double desired_error_rate)
{
	double min = 0.00015;
	double max = 6000;
	while(1) {
		double mid = (min + max) / 2;
		double error_rate = 0;
		double len;

		qscore_calibrate(mid, qval_to_qval, num_quality_scores);
		for(unsigned i = 0; i < num_quality_scores; i++)
			error_rate += qscore_ratio_distr[i] *
				      qscore_to_error_rate(qval_to_qval[i]);

		len = 0;
		if (error_rate > desired_error_rate) {
			len = max - mid;
			max = mid;
		} else if (error_rate < desired_error_rate) {
			len = mid - min;
			min = mid;
		}
		
		if (fabs(error_rate - desired_error_rate) < 0.0001 || len < 0.00001)
			return;
	}
}

// 
// Adjusts a row of the distribution matrix to take into account the new quality
// score mapping @qval_to_qval.
//
static void adjust_dist_matrix_row(uint64_t dist_matrix_row[],
				   const unsigned qval_to_qval[],
				   unsigned num_seq_bases,
				   unsigned num_quality_scores)
{
	unsigned row_len = num_seq_bases * num_quality_scores;
	uint64_t orig_dist_matrix_row[row_len];
	memcpy(orig_dist_matrix_row, dist_matrix_row, row_len * sizeof(uint64_t));
	memset(dist_matrix_row, 0, row_len * sizeof(uint64_t));

	for (unsigned old_qval = 0; old_qval < num_quality_scores; old_qval++) {
		unsigned new_qval = qval_to_qval[old_qval];
		for (unsigned seq_base = 0; seq_base < num_seq_bases; seq_base++) {
			dist_matrix_row[seq_base * num_quality_scores + new_qval] += 
				orig_dist_matrix_row[seq_base * num_quality_scores + new_qval];
		}
	}
}

// 
// Transforms a row of the distribution matrix to achieve an average error rate
// of approximately @desired_error_rate.
//
static void user_error_rate_transform(uint64_t dist_matrix_row[],
				      unsigned num_seq_bases,
				      unsigned num_quality_scores,
				      double desired_error_rate,
				      unsigned bin_ref_base)
{
	unsigned row_len = num_seq_bases * num_quality_scores;
	uint64_t qscore_distr[num_quality_scores];
	double qscore_ratio_distr[num_quality_scores];
	unsigned qval_to_qval[num_quality_scores];
	unsigned i;
	uint64_t total_correct_sum = 0;
	uint64_t total_error_sum = 0;
	uint64_t total_sum;

	ZERO_ARRAY(qscore_distr);

	for (i = 0; i < row_len; i++) {
		unsigned qscore = i % num_quality_scores;
		unsigned seq_base = i / num_quality_scores;
		uint64_t freq = dist_matrix_row[i];
		qscore_distr[qscore] += freq;
		if (seq_base == bin_ref_base)
			total_correct_sum += freq;
		else
			total_error_sum += freq;
	}
	total_sum = total_correct_sum + total_error_sum;

	for (i = 0; i < num_quality_scores; i++)
		qscore_ratio_distr[i] = double(qscore_distr[i]) / double(total_sum);
	
	compute_qval_to_qval(qval_to_qval, qscore_ratio_distr,
			     num_quality_scores, desired_error_rate);

	adjust_dist_matrix_row(dist_matrix_row, qval_to_qval, num_seq_bases,
			       num_quality_scores);
}

/* Load a line of the distribution matrix. */
void load_dist_matrix_line(char *line, Profile &_profile)
{
	BaseCallingProfile &profile = static_cast<BaseCallingProfile&>(_profile);
	char *ref_base_str;
	unsigned bin_ref_base;
	int cycle;
	uint64_t *dist_matrix_row;
	uint64_t **dist_called_base_matrix_section;
	size_t row_len;

	ref_base_str = strtok(line, "\t ");
	if (!ref_base_str)
		return;

	bin_ref_base = dna_char_to_bin(ref_base_str[0]);
	if (bin_ref_base >= 4)
		return;

	cycle = profile.get_cycle(strtok(NULL, "\t "));
	if (cycle == -1)
		return;

	dist_matrix_row = profile.dist_matrix[cycle - 1][bin_ref_base];
	row_len = profile.num_seq_bases * profile.num_quality_scores;

	if (!load_matrix_row(dist_matrix_row, row_len))
		goto invalid;


	if (profile.params.error_rate != -1.0) {
		user_error_rate_transform(dist_matrix_row,
					  profile.num_seq_bases,
					  profile.num_quality_scores,
					  profile.params.error_rate,
					  bin_ref_base);
	}

	if (profile.params.subst_error_algo == ALGO_QTRANS) {
		dist_called_base_matrix_section = profile.dist_called_base_matrix[cycle - 1][bin_ref_base];
		for (int j = 0; j < profile.num_quality_scores; j++)
			memset(dist_called_base_matrix_section[j], 0, profile.num_seq_bases * sizeof(uint64_t));
		for (int i = 0; i < profile.num_seq_bases; i++)
			for (int j = 0; j < profile.num_quality_scores; j++)
				dist_called_base_matrix_section[j][i] += dist_matrix_row[i * profile.num_quality_scores + j];

		for (int j = 0; j < profile.num_quality_scores; j++) {
			if (!prepare_uint64_probability_array(dist_called_base_matrix_section[j],
							      profile.num_seq_bases,
							      bin_ref_base))
				goto invalid;
		}
	}

	if (!prepare_uint64_probability_array(dist_matrix_row, row_len,
					      profile.num_quality_scores * (bin_ref_base + 1) - 1))
		goto invalid;
	return;
invalid:
	fatal_error("Invalid base-calling profile; please check "
		    "that it is in the right format!");
}

/* Load a line of the quality-transition matrix. */
void load_qtrans_matrix_line(char *line, Profile &_profile)
{
	BaseCallingProfile &profile = static_cast<BaseCallingProfile&>(_profile);
	int cycle;
	char *cycle_str = strtok(line, "\t ");
	int prev_qscore;
	size_t row_len;
	uint64_t *qtrans_matrix_row;
	const char *prev_qscore_str;

	cycle = profile.get_cycle(cycle_str);
	if (cycle == -1)
		return;

	prev_qscore_str = strtok(NULL, "\t ");
	if (!get_int(prev_qscore_str, prev_qscore))
		return;

	if (prev_qscore < 0 || prev_qscore >= profile.num_quality_scores)
		goto invalid;

	qtrans_matrix_row = profile.qtrans_matrix[cycle - 1][prev_qscore];
	row_len = profile.num_quality_scores;


	if (!load_matrix_row(qtrans_matrix_row, row_len))
		goto invalid;

	if (profile.params.error_rate != -1.0) {
		user_error_rate_transform(qtrans_matrix_row, 1,
					  profile.num_quality_scores,
					  profile.params.error_rate, 0);
	}

	prepare_uint64_probability_array(qtrans_matrix_row, row_len,
					 profile.num_quality_scores - 1);
	return;
invalid:
	fatal_error("Invalid base-calling profile; please check "
		    "that it is in the right format!");
}

/* 
 * Load the header of the base-calling profile; it provides some information
 * about the dimensionality of the matrices.
 */
void BaseCallingProfile::load_dimensions_header()
{
	InputStream in_file(filename);
	char *line = NULL;
	size_t n;
	while (in_file.getline(&line, &n) != -1)
		if (sscanf(line, "# Dimensions: Ref_base_number %d, Cycle_number %d, "
			   "Seq_base_number %d, Quality_number %d", 
			   &num_ref_bases, &num_cycles,
			   &num_seq_bases, &num_quality_scores) == 4)
			goto found_dimensions;

	fatal_error("Could not find '#Dimensions' line in base-calling "
		    "profile \"%s\"!\nPlease make sure the base-calling "
		    "profile is in the correct format.",
		    filename.c_str());
found_dimensions:
	if (num_ref_bases <= 0 || num_cycles <= 0 || 
	    num_seq_bases <= 0 || num_quality_scores <= 0) 
		fatal_error("At least one of num_ref_bases, num_cycles, "
			    "num_seq_bases, and/or num_quality_scores was invalid.\n"
			    "Please check the base-calling profile \"%s\"!",
			    filename.c_str());

	if (2 * params.read_len > num_cycles)
		fatal_error("The base-calling profile \"%s\" supports at most %d\n"
			    "base-pair reads.  Please lower the read length, or choose a different\n"
			    "base-calling profile.",
			    filename.c_str(), num_cycles / 2);

	info("Loaded base-calling profile header:\n");
	info("   num_ref_bases      = %d\n", num_ref_bases);
	info("   num_cycles         = %d\n", num_cycles);
	info("   num_seq_bases      = %d\n", num_seq_bases);
	info("   num_quality_scores = %d\n", num_quality_scores);
	free(line);
}

/*
 * Construct the base-calling profile based on the simulation parameters.
 */
BaseCallingProfile::BaseCallingProfile(const SimulationParameters & params)
	: Profile		 (params),	
	  num_ref_bases		 (0),
	  num_cycles		 (0),
	  num_seq_bases		 (0),
	  num_quality_scores	 (0),
	  dist_matrix		 (NULL),
	  dist_called_base_matrix(NULL),
	  qtrans_matrix		 (NULL)
{
	if (params.base_calling_profile_filename)
		filename = params.base_calling_profile_filename;
	else
		filename = DEFAULT_BASE_CALLING_PROFILE;

	info("Loading base-calling profile %s\n", filename.c_str());

	load_dimensions_header();

	dist_matrix = (uint64_t***) new_matrix(3, params.read_len * 2,
					       num_ref_bases, 
					       num_seq_bases * num_quality_scores);
	if (params.subst_error_algo == ALGO_QTRANS) {
		qtrans_matrix = (uint64_t***)new_matrix(3, params.read_len * 2,
					   	        num_quality_scores,
							num_quality_scores);
		dist_called_base_matrix = (uint64_t****)new_matrix(
						4, params.read_len * 2,
						num_ref_bases, num_quality_scores,
						num_seq_bases);
	}

	info("Loading DistMatrix of base-calling profile\n");
	for_line_in_matrix("[DistMatrix]", load_dist_matrix_line);

	if (params.subst_error_algo == ALGO_QTRANS) {
		info("Loading QTransMatrix of base-calling profile\n");
		for_line_in_matrix("[QTransMatrix]", load_qtrans_matrix_line);
	}

	info("Done loading Base-Calling profile.\n");
}

BaseCallingProfile::~BaseCallingProfile()
{
	delete_matrix(dist_matrix, 3, params.read_len * 2,
		      num_ref_bases, num_seq_bases * num_quality_scores);

	delete_matrix(dist_called_base_matrix, 4, params.read_len * 2,
		      num_ref_bases, num_quality_scores, num_seq_bases);

	delete_matrix(qtrans_matrix, 3, params.read_len * 2,
		      num_quality_scores, num_quality_scores);
}

/*
 * Call a base, based on the base-calling profile and the specified parameters.
 *
 * @cycle: The current cycle (1 <= cycle <= 2 * params.read_len).
 * @raw_base:  The "actual" base.  This usually will be the same as in the
 * 		reference sequence, but this is actually after the indel error 
 * 		step, so it could possibly just be a random base from an insert.
 * @prev_qscore: The quality score in the previous cycle.
 * @qscore:	Reference to an int into which the called quality score will be
 * 		returned.
 *
 * The return value is the called base.
 */
char BaseCallingProfile::call(int cycle, char raw_base, int prev_qscore,
			      int &qscore, Random &rgen) const
{
	int bin_raw_base = dna_char_to_bin(raw_base);
	char called_base;
	size_t location;
	const uint64_t *dist_matrix_row = dist_matrix[cycle - 1][bin_raw_base];
	if (params.subst_error_algo == ALGO_DIST || cycle == 1 || cycle == params.read_len + 1) {
		location = rgen.search_location(dist_matrix_row,
						num_quality_scores * num_seq_bases);
		called_base = dna_bin_to_char(location / num_quality_scores);
		qscore = location % num_quality_scores;
	} else {
		unsigned bin_called_base;
		const uint64_t *qtrans_matrix_row;
		const uint64_t *dist_called_base_matrix_row;

		assert(params.subst_error_algo == ALGO_QTRANS);

		qtrans_matrix_row = qtrans_matrix[cycle - 1][prev_qscore];
		qscore = rgen.search_location(qtrans_matrix_row,
					      num_quality_scores);
		dist_called_base_matrix_row = dist_called_base_matrix[cycle - 1][bin_raw_base][qscore];
		bin_called_base = rgen.search_location(dist_called_base_matrix_row,
						       num_seq_bases);
		called_base = dna_bin_to_char(bin_called_base);
	}
	return called_base;
}
