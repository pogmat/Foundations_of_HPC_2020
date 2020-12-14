/*
 * pgm_bin.h
 * Simple header to read and write .pgm binary files (magic value "P5").
 *
 * Author: Matteo Poggi <matteo.poggi.fi@gmail.com>
 */ 

#ifndef _PGM_BIN_H
#define _PGM_BIN_H

#include <stdint.h>

typedef uint8_t byte;
typedef uint16_t dbyte;

typedef struct {
	char magic[3];
	unsigned int width;
	unsigned int height;
	unsigned int maximum_value;
	byte* data;
} pgm_file;

int read_pgm(const char* const, pgm_file* const);
int write_pgm(const char* const, const pgm_file* const);
void free_pgm(pgm_file* const);

#endif
