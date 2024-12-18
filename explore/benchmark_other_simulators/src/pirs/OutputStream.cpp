#include "OutputStream.h"
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/* See OutputStream.h for documentation. */

enum OutputType OutputStream::default_output_type = TEXT;

enum OutputType OutputStream::str_to_output_type(const char *str)
{
	if (strcmp(str, "txt") == 0 || strcmp(str, "text") == 0)
		return TEXT;
	else if (strcmp(str, "gz") == 0 || strcmp(str, "gzip") == 0)
		return GZIP;
	else
		fatal_error("-c (--output-file-type) option must be "
			    "given the argument \"text\" or \"gzip\"!");
}

const char *OutputStream::get_default_file_type_str()
{
	switch (default_output_type) {
	case TEXT:
		return "text";
	case GZIP:
		return "gzip";
	default:
		assert(0);
	}
}

const char *OutputStream::get_file_type_suffix()
{
	switch (s_type) {
	case TEXT:
		return "";
	case GZIP:
		return ".gz";
	default:
		assert(0);
	}
}

void OutputStream::open(const char *filename, enum OutputType type)
{
	bool is_stdout = (strcmp(filename, "-") == 0);
	free(s_filename);
	
	if (type == DEFAULT)
		s_type = default_output_type;
	else
		s_type = type;

	if (strcmp(filename, "-") == 0) {
		s_filename = xstrdup("(stdout)");
		is_stdout = true;
	} else {
		s_filename = xstrdup2(filename, get_file_type_suffix());
		is_stdout = false;
	}

	//info("Opening output file \"%s\"\n", s_filename);

	switch (s_type) {
	case TEXT:
		if (is_stdout)
			text.s_fp = stdout;
		else
			text.s_fp = fopen(s_filename, "wb");
		if (!text.s_fp)
			goto fail;
		break;
	case GZIP:
		if (is_stdout)
			gzip.s_fp = gzdopen(STDOUT_FILENO, "wb");
		else
			gzip.s_fp = gzopen(s_filename, "wb");
		if (!gzip.s_fp)
			goto fail;
		break;
	default:
		assert(0);
	}
	return;
fail:
	fatal_error_with_errno("Failed to open the file \"%s\" for writing",
			       s_filename);
}


void OutputStream::write(const void *buf, size_t len)
{
	switch (s_type) {
	case TEXT:
		if (fwrite(buf, 1, len, text.s_fp) != len)
			goto fail;
		break;
	case GZIP:
		if (gzwrite(gzip.s_fp, buf, len) != (int)len)
			goto fail;
		break;
	default:
		assert(0);
	}
	return;
fail:
	fatal_error_with_errno("Error writing to the file \"%s\"",
			       s_filename);
}

void OutputStream::write_header()
{
	this->printf("# File \"%s\": generated at %s by the command:\n"
		     "#     %s\n",
		     s_filename, get_timestamp(), get_command_line());
}

void OutputStream::vprintf(const char *format, va_list va)
{
	switch (s_type) {
	case TEXT:
		if (vfprintf(text.s_fp, format, va) < 0)
			goto fail;
		break;
	case GZIP: {
			int len;
			va_list va_save;
			if (!gzip.s_aux_buf) {
				gzip.s_aux_buf_sz = 128;
				gzip.s_aux_buf = (char*)xmalloc(gzip.s_aux_buf_sz);
			}
			while (1) {
				va_copy(va_save, va);
				len = vsnprintf(gzip.s_aux_buf,
						gzip.s_aux_buf_sz,
						format, va);
				va_end(va);
				if (len < gzip.s_aux_buf_sz)
					break;
				va_copy(va, va_save);
				va_end(va_save);
				gzip.s_aux_buf_sz *= 2;
				gzip.s_aux_buf = (char*)xrealloc(gzip.s_aux_buf,
							         gzip.s_aux_buf_sz);
			}
			if (gzwrite(gzip.s_fp, gzip.s_aux_buf, len) != len) {
				int errnum;
				const char *err_str = gzerror(gzip.s_fp, &errnum);
				if (errnum == Z_ERRNO)
					goto fail;
				fatal_error("zlib error writing to the file "
					    "\"%s\": %s", s_filename, err_str);
			}
		}
		break;
	default:
		assert(0);
	}
	return;
fail:
	fatal_error_with_errno("Error writing to the file \"%s\"",
			       s_filename);
}

void OutputStream::printf(const char *format, ...)
{
	va_list va;
	va_start(va, format);
	this->vprintf(format, va);
	va_end(va);
}

void OutputStream::putc(int c)
{
	switch (s_type) {
	case TEXT:
		if (::putc(c, text.s_fp) == EOF)
			goto fail;

		break;
	case GZIP:
		if (::gzputc(gzip.s_fp, c) == -1)
			goto fail;
		break;
	default:
		assert(0);
	}
	return;
fail:
	fatal_error_with_errno("Error writing to the file \"%s\"",
			       s_filename);
}

void OutputStream::close()
{
	if (!is_open())
		return;

	//info("Closing output file \"%s\"\n", s_filename);
	switch (s_type) {
	case TEXT:
		if (text.s_fp && fclose(text.s_fp) != 0)
			goto fail;
		break;
	case GZIP: {
			int errnum;
			const char *err_str;
			if (gzip.s_fp && (errnum = gzclose(gzip.s_fp)) != Z_OK) {
				if (errnum == Z_ERRNO)
					goto fail;
				err_str = gzerror(gzip.s_fp, &errnum);
				fatal_error("zlib error closing the file "
					    "\"%s\": %s", s_filename, err_str);

			}
		}
		free(gzip.s_aux_buf);
		gzip.s_aux_buf = NULL;
		break;
	default:
		assert(0);
	}
	free(s_filename);
	s_filename = NULL;
	return;
fail:
	fatal_error_with_errno("Error closing the file \"%s\"", s_filename);
}

