#include "IndelProfile.h"
#include "SimulationParameters.h"
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using std::string;

//
// Load a line of the InDel profile.
//
// Each line of the matrix contains the data for a cycle.  There is one number
// for each possible indel length (-profile.max_indel_len to
// +profile.max_indel_len), other than 0, which is omitted.  The 0 entry
// represents no indel and can be determined by taking the total number of reads
// and subtracting the number of indels from the total number of reads.
//
static void load_indel_error_matrix_line(char *line, Profile &_profile)
{
	//[InDel]
	//Cycle   -3      -2      -1      1       2       3
	//1       0       6       2       4       0       0
	//2       0       2       13      32      5       0
	//3       0       5       224     31      3       0
	//4       0       1       278     34      7       0
	//.....
	//<<END
	
	IndelProfile &profile = static_cast<IndelProfile&>(_profile);
	const SimulationParameters &params = profile.params;
	char *tmp;
	char *cycle_str = strtok(line, "\t ");

	if (!cycle_str)
		return;

	int cycle = strtol(cycle_str, &tmp, 10);
	if (tmp == cycle_str)
		return;
    		
	if (cycle > profile.num_cycles / 2 && 
			cycle <= profile.num_cycles + params.read_len)
		cycle = params.read_len + (cycle - profile.num_cycles / 2);
	else if (cycle < 1 || cycle > params.read_len)
		return;

	uint64_t *indel_matrix_row = profile.indel_error_matrix[cycle - 1];
	uint64_t *row_center = &indel_matrix_row[profile.max_indel_len];
	uint64_t indel_sum = 0;

	for (int i = -profile.max_indel_len; i <= profile.max_indel_len; i++) {

		if (i == 0)
			continue;

		char *n_str = strtok(NULL, "\t ");
		if (!n_str)
			goto invalid;
		uint64_t n = strtoull(n_str, &tmp, 10);
		if (tmp == n_str)
			goto invalid;
		if (i > 0) {
			// insertion
			if (cycle <= params.read_len)
				profile.read_1_ins_total_len += n * i;
			else
				profile.read_2_ins_total_len += n * i;
		} else {
			// deletion
			if (cycle <= params.read_len)
				profile.read_1_del_total_len += n * -i;
			else
				profile.read_2_del_total_len += n * -i;
		}

		row_center[i] = n;
		indel_sum += n;
	}
	if (strtok(NULL, " \t"))
		goto invalid;

	if (cycle < params.read_len)
		*row_center = profile.read_1_count - indel_sum;
	else
		*row_center = profile.read_2_count - indel_sum;

	if (!prepare_uint64_probability_array(indel_matrix_row,
					 2 * profile.max_indel_len + 1,
					 profile.max_indel_len))
		goto invalid;
	return;
invalid:
	fatal_error("Invalid InDel error profile; please check "
		    "that it is in the right format!");
}

// Load a line from the header of the indel error profile
static void load_indel_info_line(char *line, Profile &_profile)
{
	IndelProfile &profile = static_cast<IndelProfile&>(_profile);
	const SimulationParameters &params = profile.params;
	const char *filename = profile.filename.c_str();
	uint64_t n;
	int n_int;
	
	if (sscanf(line, "Read_Length = %d", &n_int) == 1) {
		profile.num_cycles = n_int * 2;
		if (params.read_len > n_int) {
			fatal_error("The selected InDel error profile "
				    "(\"%s\") can only handle up to "
				    "%d bp reads.\nPlease choose a "
				    "lower read length, or a different "
				    "InDel error profile.\n", 
				    filename, n_int);
		}
	} else if (sscanf(line, "Read_1_Count = %lu", &n) == 1)
		profile.read_1_count = n;
	else if (sscanf(line, "Read_2_Count = %lu", &n) == 1)
		profile.read_2_count = n;
	else if (sscanf(line, "MaxInDel_Length = %d", &n_int) == 1)
		profile.max_indel_len = n_int;
}

/* Construct the indel profile from the simulation parameters. */
IndelProfile::IndelProfile(const SimulationParameters &params)
	: Profile(params),
	  indel_error_matrix	(NULL),
	  max_indel_len		(0),
	  num_cycles 		(0),
	  read_1_count		(0),
	  read_2_count		(0),
	  read_1_ins_total_len	(0),
	  read_1_del_total_len	(0),
	  read_2_ins_total_len 	(0),
	  read_2_del_total_len	(0)

{
	if (params.indel_profile_filename)
		filename = params.indel_profile_filename;
	else
		filename = DEFAULT_INDEL_PROFILE;
	
	info("Loading InDel error profile %s\n", filename.c_str());
	
	for_line_in_matrix("[Info]", load_indel_info_line);

	if (num_cycles == 0 || max_indel_len == 0) {
		fatal_error("The InDel profile \"%s\"\nis invalid!  Without "
			    "an appropriate InDel profile, indel errors cannot "
			    "be simulated.\nUse the --no-indels option "
			    "if you do not want to simulate indels,\nor make "
			    "sure your InDel-error profile is in the right "
			    "format.", filename.c_str());
	}

	indel_error_matrix = (uint64_t**)new_matrix(2, num_cycles, 
						  max_indel_len * 2 + 1);

	for_line_in_matrix("[InDel]", load_indel_error_matrix_line);

	info("Indel error profile statistics:\n");
	info("    max_indel_len = %d\n", max_indel_len);
	info("    num_cycles    = %d\n", num_cycles);
	info("    read_1_count  = %lu\n", read_1_count);
	info("    read_2_count  = %lu\n", read_2_count);

	double read_1_ins_rate = (double)read_1_ins_total_len /
				 (double)(read_1_count * params.read_len);

	double read_1_del_rate = (double)read_1_del_total_len /
				 (double)(read_1_count * params.read_len);

	double read_2_ins_rate = (double)read_2_ins_total_len /
				 (double)(read_2_count * params.read_len);

	double read_2_del_rate = (double)read_2_del_total_len /
				 (double)(read_2_count * params.read_len);

	info("    Insertion-bp rate of %d-bp read 1 = %.6f%%\n",
	     params.read_len, read_1_ins_rate * 100.0);

	info("    Deletion-bp rate of %d-bp read 1  = %.6f%%\n",
	     params.read_len, read_1_del_rate * 100.0);

	info("    Insertion-bp rate of %d-bp read 2 = %.6f%%\n",
	     params.read_len, read_2_ins_rate * 100.0);

	info("    Deletion-bp rate of %d-bp read 2  = %.6f%%\n",
	     params.read_len, read_2_del_rate * 100.0);

	info("Done loading InDel error profile.\n");
}
