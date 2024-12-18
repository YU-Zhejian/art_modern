#ifndef _INDEL_PROFILE_H
#define _INDEL_PROFILE_H

#include "Profile.h"
#include "Random.h"

#include <inttypes.h>

/*
 * Profile for simulating indels in the reads.
 */
class IndelProfile : public Profile {
public:
	uint64_t  **indel_error_matrix;
	int	    max_indel_len;
	int	    num_cycles;
	uint64_t    read_1_count;
	uint64_t    read_2_count;
	uint64_t    read_1_ins_total_len;
	uint64_t    read_1_del_total_len;
	uint64_t    read_2_ins_total_len;
	uint64_t    read_2_del_total_len;

	/* Construct the indel profile from the simulation parameters. */
	IndelProfile(const SimulationParameters &params);

	~IndelProfile() {
		delete_matrix(indel_error_matrix, 2, num_cycles,
			      2 * max_indel_len + 1);
	}

	/*
	 * Get the length of an indel, given a cycle number.  0 means no indel;
	 * negative means a deletion, and positive means an insertion.
	 */
	int get_indel_len(int cycle, Random &rgen) const {
		size_t location = rgen.search_location(indel_error_matrix[cycle],
						       max_indel_len * 2 + 1);
		return location - max_indel_len;
	}
};

#define DEFAULT_INDEL_PROFILE	     PKGDATADIR"/InDel_Profiles/phixv2.InDel.matrix"

#endif /* ifndef _INDEL_PROFILE_H */
