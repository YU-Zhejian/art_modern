#ifndef _PROFILE_H
#define _PROFILE_H


#include <string>

class SimulationParameters;
class Profile;

using std::string;

typedef void (*matrix_processor_func)(char *line, Profile &profile);

/*
 * Functionality shared by all profiles (the filename of the profile, and the
 * parameters from which the profile was constructed)
 */
class Profile {
public:
	string filename;
	const SimulationParameters &params;

	Profile(const SimulationParameters &__params) 
		: params(__params)
	{ }
protected:
	void for_line_in_matrix(const char *matrix_tag, 
				matrix_processor_func f);
};


#endif /* ifndef _PROFILE_H */
