#include "Profile.h"
#include "InputStream.h"
#include <string.h>
#include <stdlib.h>

/*
 * Call @f on each line in the Profile's file that is within the matrix tagged
 * by @matrix_tag, such as "[DistMatrix]".
 */
void Profile::for_line_in_matrix(const char *matrix_tag,
				 matrix_processor_func f)
{

	InputStream in_file(filename);
	size_t tag_len = strlen(matrix_tag);
	char tag[tag_len + 2];
	strcpy(tag, matrix_tag);
	strcat(tag, "\n");
	bool inside_matrix = false;
	char *line = NULL;
	size_t n;
	while (in_file.getline(&line, &n) != -1) {
		if (strcmp(line, tag) == 0) {
			inside_matrix = true;
		} else {
			if (inside_matrix) {
				if (!*line || line[0] == '#')
					continue;
				if (strcmp(line, "<<END\n") == 0)
					break;
				f(line, *this);
			}
		}
	}
	free(line);
}
