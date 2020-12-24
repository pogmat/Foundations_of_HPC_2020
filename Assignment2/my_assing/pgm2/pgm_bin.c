/*
 * pgm_bin.c
 * Simple source to read and write .pgm binary files (magic value "P5").
 *
 * Author: Matteo Poggi <matteo.poggi.fi@gmail.com>
 */ 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pgm_bin.h"

#if ((0x100 & 0xf) == 0x0)
#define ENDIAN_SWAP 1
#else
#define ENDIAN_SWAP 0
#endif

static void skip_comment(FILE * const fp)
{
	int ch;
	int inside_comment = 0;
	
	while ((ch = fgetc(fp)) != EOF) {
		switch (ch) {
		case ' ':
		case '\t':
			break;
		case '#':
			inside_comment = 1;
			break;
		case '\n':
			inside_comment = 0;
			break;
		default:
			if (!inside_comment) {
				fseek(fp, -1, SEEK_CUR);
				return;
			}
			break;
		}
	}
}

pgm_file_t* open_pgm(const char* const filename, const char mode)
{
        pgm_file_t* pgm = (pgm_file_t*)malloc(sizeof(pgm_file_t));
	if (!pgm) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}

	pgm->image.data = NULL;
	
	switch (mode) {
	case 'w':
		pgm->file = fopen(filename, "wb");
		if (!pgm->file) {
			fprintf(stderr, "No such file: %s.\n", filename);
			return NULL;
		}
		break;
	case 'r':
		pgm->file = fopen(filename, "rb");
		if (!pgm->file) {
			fprintf(stderr, "No such file: %s.\n", filename);
			return NULL;
		}

		skip_comment(pgm->file);
		fgets(pgm->image.magic, sizeof(pgm->image.magic), pgm->file);
		if (strcmp(pgm->image.magic, "P5")) {
			fprintf(stderr, "Invalid format: only P5 format allowed.\n");
			free(pgm);
			fclose(pgm->file);
			return NULL;
		}

		skip_comment(pgm->file);
		if (!fscanf(pgm->file, "%u", &pgm->image.width)) {
			fprintf(stderr, "Invalid format: no width.\n");
			free(pgm);
			fclose(pgm->file);
			return NULL;
		}		

		skip_comment(pgm->file);
		if (!fscanf(pgm->file, "%u", &pgm->image.height)) {
			fprintf(stderr, "Invalid format: no height.\n");
			free(pgm);
			fclose(pgm->file);
			return NULL;
		}

		skip_comment(pgm->file);
		if (!fscanf(pgm->file, "%u", &pgm->image.maximum_value)) {
			fprintf(stderr, "Invalid format: no maximum value.\n");
			free(pgm);
			fclose(pgm->file);
			return NULL;
		}
		skip_comment(pgm->file);

		size_t img_size = (pgm->image.width * pgm->image.height)
			* (1 + (pgm->image.maximum_value > UINT8_MAX));
		pgm->image.data = (byte*)malloc(img_size);
		if (!pgm->image.data) {
			fprintf(stderr, "Bad allocation.\n");
			free(pgm);
			fclose(pgm->file);
			return NULL;	
		}
		break;
	default:
		fprintf(stderr, "Unknown option %c.\n", mode);
		return NULL;
		break;
	}
	
	return pgm;
}

int read_pgm(const pgm_file_t* const pgm)
{
	int two_bytes = pgm->image.maximum_value > UINT8_MAX;
	size_t img_elements = pgm->image.width * pgm->image.height;
	byte word[2];

	if (two_bytes && ENDIAN_SWAP) {
		for (size_t read_e = 0; read_e < img_elements; ++read_e) {
			if (!fread(word, sizeof(dbyte), 1, pgm->file)) {
				fprintf(stderr, "Invalid format: the file is corrupted.\n");
				return -1;
			}

			pgm->image.data[2 * read_e] = word[1];
			pgm->image.data[2 * read_e + 1] = word[0];
		}
	} else {
		if (!fread(pgm->image.data, sizeof(byte), img_elements, pgm->file)) {
			fprintf(stderr, "Invalid format: the file is corrupted.\n");
			return -1;
		}
	}

	return 0;
}

int write_pgm(const pgm_file_t* const pgm)
{
	if(!fprintf(pgm->file, "%2s\n%d %d\n%d\n",
		    pgm->image.magic,
		    pgm->image.width,
		    pgm->image.height,
		    pgm->image.maximum_value)) {
		fprintf(stderr, "Impossible to write file header.\n");
		return -1;
	}

	int two_bytes = pgm->image.maximum_value > UINT8_MAX;
	size_t img_elements = pgm->image.width * pgm->image.height;
	byte word[2];

	if (two_bytes && ENDIAN_SWAP) {
		for (size_t written_e = 0; written_e < img_elements; ++written_e) {
			//buffer = ENDIAN_ADJUST(((dbyte*)pgm->data)[written_b]);
			word[0] = pgm->image.data[2 * written_e + 1];
			word[1] = pgm->image.data[2 * written_e];
			if (!fwrite(word, sizeof(dbyte), 1, pgm->file)) {
				fprintf(stderr, "Impossible to write file data.\n");
				return -1;
			}
			
		}
	} else {
		if (!fwrite(pgm->image.data, sizeof(byte), img_elements, pgm->file)) {
			fprintf(stderr, "Impossible to write file data.\n");
			return -1;
		}
	}

	return 0;	
}

void close_pgm(pgm_file_t* const pgm)
{
	if (!pgm)
		return;

	if (pgm->image.data)
		free(pgm->image.data);

	if (pgm->file)
		fclose(pgm->file);

	free(pgm);
}
