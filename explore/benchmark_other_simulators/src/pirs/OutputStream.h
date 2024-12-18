#ifndef _OUTPUT_STREAM_H
#define _OUTPUT_STREAM_H

#include <zlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <string>
#include <stdarg.h>
#include "util.h"

using std::string;

#ifdef putc
#undef putc
#endif

/*
 * Custom output stream class to handle some of the annoying things that often
 * come up when working with files.
 *
 * open(), printf(), vprintf(), putc(), puts(), and write() methods are
 * provided.  Everything will exit the program with a (hopefully useful) error
 * message if there is an error, so the rest of the code need not worry about
 * this.
 *
 * Transparent compression is supported. Specify the output type in the open(),
 * or set the default type with  set_default_output_type().  The default is
 * uncompressed.
 *
 * The destructor will automatically close the file, with error checking.
 */
class OutputStream {
private:
	static enum OutputType default_output_type;

	enum OutputType s_type;
	union {
		struct {
			gzFile s_fp;
			char *s_aux_buf;
			int s_aux_buf_sz;
		} gzip;
		struct {
			FILE *s_fp;
		} text;
	};
	const char *get_file_type_suffix();
	static enum OutputType str_to_output_type(const char *str);
public:
	char *s_filename;

	OutputStream() { 
		memset(this, 0, sizeof(*this));
	}

	OutputStream(const char *filename, enum OutputType type = DEFAULT) {
		memset(this, 0, sizeof(*this));
		open(filename, type);
	}

	~OutputStream() {
		this->close();
	}

	void open(const char *filename, enum OutputType type = DEFAULT);
	void open(const string &filename, enum OutputType type = DEFAULT) {
		open(filename.c_str(), type);
	}
	bool is_open() {
		return s_filename != NULL;
	}
	void write(const void *buf, size_t nbyte);
	void write_header();
	void printf(const char *format, ...) __format(printf, 2, 3);
	void vprintf(const char *format, va_list va) __format(printf, 2, 0);
	void puts(const char *str) { 
		this->write(str, strlen(str));
	}
	void puts(const string &str) {
		this->write(str.c_str(), str.length());
	}
	void putc(int c);
	void close();

	static void set_default_output_type(const char *str) {
		default_output_type = str_to_output_type(str);
	}
	static void set_default_output_type(enum OutputType type) {
		default_output_type = type;
	}
	static bool default_output_type_is_compressed() {
		return default_output_type != TEXT;
	}
	static const char *get_default_file_type_str();
};

#endif /* _OUTPUT_STREAM_H */
