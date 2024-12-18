#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <string>

#include "util.h"
#include "pirs.h"
using std::string;

void program_info()
{
	const char *info_str = 
"Program:       " PROGRAM "\n"
"Version:       " VERSION "\n"
"Author:        " AUTHOR "\n"
"Contact:       " CONTACT "\n"
"Compile Date:  " __DATE__ " time: " __TIME__ "\n"
	;
	fputs(info_str, stdout);
}

static void pirs_usage()
{
	const char *usage_str = 
"Usage: pirs <command> [option]\n"
"    diploid     generate diploid genome.\n"
"    simulate    simulate Illumina reads.\n"
	;
	fputs(usage_str, stdout);
}

const char *timestamp;
const char *command_line;

void set_timestamp()
{
	static string timestamp_string;

	time_t now = time(NULL);
	char *p = ctime(&now);
	p[strlen(p) - 1] = '\0';
	timestamp_string = p;
	timestamp = timestamp_string.c_str();
}

void set_command_line(int argc, const char * const * argv)
{
	static string command_line_string;
	for (int i = 0; i < argc; i++) {
		command_line_string += argv[i];
		if (i != argc - 1)
			command_line_string += ' ';
	}
	command_line = command_line_string.c_str();
}


int main(int argc, char *argv[])
{
	set_timestamp();
	set_command_line(argc, argv);

	if (argc < 2) {
		program_info();
		putchar('\n');
		pirs_usage();
		return 2;
	}
#ifdef ENABLE_PIRS_DIPLOID
	if (strcmp("diploid", argv[1]) == 0) {
		pirs_diploid(argc, argv);
		return 0;
	}
#endif
#ifdef ENABLE_PIRS_SIMULATE
	if (strcmp("simulate",argv[1]) == 0) {
		pirs_simulate(argc, argv);
		return 0;
	}
#endif

	program_info();
	putchar('\n');
	pirs_usage();
	return 2;
}

