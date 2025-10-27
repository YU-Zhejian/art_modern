/**
 * @file getdelim.h
 * @brief NetBSD getdelim and getline declarations.
 */

#ifndef AM_NETBSD_GETDELIM_H
#define AM_NETBSD_GETDELIM_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

ssize_t getdelim(char** buf, size_t* bufsiz, int delimiter, FILE* fp);
ssize_t getline(char** buf, size_t* bufsiz, FILE* fp);

#if defined(__cplusplus)
}
#endif
#endif /* AM_NETBSD_GETDELIM_H */
