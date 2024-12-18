#ifndef _GC_BIAS_PROFILE_H
#define _GC_BIAS_PROFILE_H

#include "Profile.h"
#include "SimulationParameters.h"
class Random;

/*
 * Profile for simulating GC bias in the reads. */
class GCBiasProfile : public Profile {
private:
	// (2**64 - 1) is 100% abundance.
	uint64_t gc_bias_abundance[101];
	unsigned char char_to_gc_count_tab[256];

public:
	GCBiasProfile(const SimulationParameters &params);
	bool accept_insert(const char *insert, size_t insert_len,
			   Random &rgen) const;
};

#define DEFAULT_GC_BIAS_PROFILE_100  PKGDATADIR"/GC-depth_Profiles/humNew.gcdep_100.dat"
#define DEFAULT_GC_BIAS_PROFILE_150  PKGDATADIR"/GC-depth_Profiles/humNew.gcdep_150.dat"
#define DEFAULT_GC_BIAS_PROFILE_200  PKGDATADIR"/GC-depth_Profiles/humNew.gcdep_200.dat"

#endif /* ifndef _GC_BIAS_PROFILE_H */
