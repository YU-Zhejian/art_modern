#include <inttypes.h>
#include <string>
#include <vector>

#include "GCBiasProfile.h"
#include "Random.h"

using std::vector;

// Load the GC bias profile from the file specified in the simulation
// parameters, or the default GC bias profile if none was specified.
GCBiasProfile::GCBiasProfile(const SimulationParameters &params)
	: Profile(params)
{
	ZERO_ARRAY(char_to_gc_count_tab);
	char_to_gc_count_tab[(unsigned char)'G'] = 1;
	char_to_gc_count_tab[(unsigned char)'C'] = 1;
	char_to_gc_count_tab[(unsigned char)'g'] = 1;
	char_to_gc_count_tab[(unsigned char)'c'] = 1;
	if (params.gc_bias_profile_filename) {
		filename = params.gc_bias_profile_filename;
	} else {
		// The default GC bias profile is determined based on the insert
		// length.
		int window_size = params.insert_len_mean;

		if (window_size < 125)
			filename = DEFAULT_GC_BIAS_PROFILE_100;
		else if (window_size < 175)
			filename = DEFAULT_GC_BIAS_PROFILE_150;
		else
			filename = DEFAULT_GC_BIAS_PROFILE_200;
	}

	info("Loading GC bias profile %s\n", filename.c_str());
	
	InputStream in_file(filename);

	char *line = NULL;
	size_t n;
	vector<double> gc_ratio_vec;
	vector<double> depth_vec;
	while (in_file.getline(&line, &n) != -1) {
		if (*line == '\0' || *line == '#' || *line == ' ')
			continue;

		double gc_ratio;
		double relative_coverage;

		//#GC%    RefCnt  DepthCnt       Mean	SmoothedMean    ...
		//0.5     2285822 1371         0.334179  0.334179   ....
		//............
		if (sscanf(line, "%lg %*s %*s %*s %lg", &gc_ratio,
			   &relative_coverage) != 2)  {
			fatal_error("The GC bias profile is in an invalid "
				    "format!  Use --no-gc-bias if you do\n"
				    "not want to simulate GC content bias. "
				    "Otherwise, make sure your GC bias profile\n"
				    "is in the right format.");
		}

		gc_ratio_vec.push_back(gc_ratio);
		depth_vec.push_back(relative_coverage);
	}
	free(line);
	
	//find max depth
	double max_depth = max_array(&depth_vec[0], depth_vec.size());

	double gc_bias_abundance_ratio[101];
	
	//convert to GC abundance
	for(size_t i = 1; i < gc_ratio_vec.size(); i++)
		for(int j = int(gc_ratio_vec[i - 1]); j < int(gc_ratio_vec[i]); j++)
			gc_bias_abundance_ratio[j] = depth_vec[i - 1] / max_depth;

	for(int i = (int)gc_ratio_vec[gc_ratio_vec.size() - 1]; i <= 100; i++)
		gc_bias_abundance_ratio[i] = depth_vec[depth_vec.size() - 1] / max_depth;

	info("GC bias profile (GC content, relative abundance), "
	     "showing every 10%%\n");
	for (int i = 0; i < 101; i++) {
		if (i % 10 == 0)
			info("    %d%%\t%g\n", i, gc_bias_abundance_ratio[i]);
		gc_bias_abundance[i] = ~0ULL * gc_bias_abundance_ratio[i];
	}

	info("Done loading GC bias profile.\n");
}

// Accepts (returns %true) or rejects (returns %false) an insert based on its GC
// contentt, and on the GC bias profile.
bool GCBiasProfile::accept_insert(const char *insert, size_t insert_len,
				  Random &rgen) const
{
	size_t gc_count = 0;
	for (size_t i = 0; i < insert_len; i++)
		gc_count += char_to_gc_count_tab[(unsigned char)insert[i]];
	
	size_t gc_percent = (100 * gc_count) / insert_len;

	// decide whether ignore this insert based on the GC content.
	return (rgen.next_uint64() <= gc_bias_abundance[gc_percent]);
}
