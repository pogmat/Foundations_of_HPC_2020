/*
 * pgm_bin.h
 * Simple header to read and write .pgm binary files (magic value "P5").
 *
 * Author: Matteo Poggi <matteo.poggi.fi@gmail.com>
 */ 

#ifndef _PGM_BIN_H
#define _PGM_BIN_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

typedef uint8_t byte;
typedef uint16_t dbyte;

typedef struct {
	char magic[3];
	unsigned int width;
	unsigned int height;
	unsigned int maximum_value;
	byte* data;
} pgm_image_t;

typedef struct {
	FILE* file;
	pgm_image_t image;
} pgm_file_t;

pgm_file_t* open_pgm(const char* const, const char);
int read_pgm(const pgm_file_t* const);
int write_pgm(const pgm_file_t* const);
void close_pgm(pgm_file_t* const);

static inline void touch_zero(const pgm_image_t* const img, size_t i) {img->data[i] = 0;}

#endif
