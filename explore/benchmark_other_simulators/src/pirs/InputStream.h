#ifndef _INPUT_STREAM_H
#define _INPUT_STREAM_H

#include <zlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <string>
#include "util.h"

using std::string;

/*
 * Custom input stream class to handle some of the annoying things that often
 * come up when working with files.
 *
 * open(), getline(), getc(), and ungetc() methods are provided.  getline()
 * works like the C library getline(), not the C++ getline().
 *
 * Transparent compression is supported; the compression type of the input
 * file is automatically determined.
 *
 * The destructor will automatically close the file.
 */
class InputStream {
private:
	enum OutputType s_type;
	union {
		struct {
			FILE *s_fp;
		} text;
		struct {
			gzFile s_fp;
		} gzip;
	};

	int gzip_status();
	int text_status();
	static enum OutputType file_type_probe(const char *filename);
public:
	char *s_filename;

	void open(const char *filename);

	void open(const string &string) { 
		open(string.c_str());
	}

	InputStream() {
		memset(this, 0, sizeof(*this));
	}

	InputStream(const char *filename) {
		memset(this, 0, sizeof(*this));
		open(filename);
	}
	InputStream(const string &filename) {
		open(filename.c_str());
	}

	bool is_open() {
		return s_filename != NULL;
	}

	ssize_t getline(char **lineptr, size_t *n);
	int getc();
	void ungetc(char c);
	void close();
	~InputStream();
};

#endif /* _INPUT_STREAM_H */
