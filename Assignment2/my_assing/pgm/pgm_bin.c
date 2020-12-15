/*
 * pgm_bin.c
 * Simple source to read and write .pgm binary files (magic value "P5").
 *
 * Author: Matteo Poggi <matteo.poggi.fi@gmail.com>
 */ 

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "pgm_bin.h"

#if ((0x100 & 0xf) == 0x0)
#define ENDIAN_SWAP 1
#else
#define ENDIAN_SWAP 0
#endif

void skip_comment(FILE * const fp)
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

int read_pgm(const char* const filename, pgm_file* const pgm)
{
	FILE* fp = fopen(filename, "rb");
	
	if (fp == NULL) {
		fprintf(stderr, "No such file: %s\n", filename);
		return -1;
	}

	skip_comment(fp);
	fgets(pgm->magic, sizeof(pgm->magic), fp);
	if (strcmp(pgm->magic, "P5")) {
		fprintf(stderr, "Invalid format: only P5 format allowed.\n");
		fclose(fp);
		return -1;
	}
	skip_comment(fp);
	if (!fscanf(fp, "%u", &pgm->width)) {
		fprintf(stderr, "Invalid format: no width.\n");
		fclose(fp);
		return -1;
	}
	skip_comment(fp);
	if (!fscanf(fp, "%u", &pgm->height)) {
		fprintf(stderr, "Invalid format: no height.\n");
		fclose(fp);
		return -1;
	}
	skip_comment(fp);
	if (!fscanf(fp, "%u", &pgm->maximum_value)) {
		fprintf(stderr, "Invalid format: no maximum value.\n");
		fclose(fp);
		return -1;
	}
	skip_comment(fp);

	int two_bytes = pgm->maximum_value > UINT8_MAX;
	size_t img_elements = pgm->width * pgm->height;
	byte word[2];

	if (two_bytes && ENDIAN_SWAP) {
		pgm->data = (byte*)malloc(img_elements * sizeof(dbyte));
		if (!pgm->data) {
			fprintf(stderr, "Bad allocation.\n");
			fclose(fp);
			return -1;
		}

		for (size_t read_e = 0; read_e < img_elements; ++read_e) {
			if (!fread(word, sizeof(dbyte), 1, fp)) {
				fprintf(stderr, "Invalid format: the file is corrupted.\n");
				fclose(fp);
				return -1;
			}

			pgm->data[2 * read_e] = word[1];
			pgm->data[2 * read_e + 1] = word[0];
		}
	} else {
		pgm->data = (byte*)malloc(img_elements * sizeof(byte));
		if (!pgm->data) {
			fprintf(stderr, "Bad allocation.\n");
			fclose(fp);
			return -1;
		}

		if (!fread(pgm->data, sizeof(byte), img_elements, fp)) {
			fprintf(stderr, "Invalid format: the file is corrupted.\n");
			fclose(fp);
			return -1;
		}
	}
	
	fclose(fp);
	
	return 0;
}

int write_pgm(const char* const filename, const pgm_file* const pgm)
{
	FILE* fp = fopen(filename, "wb");

	if (fp == NULL) {
		fprintf(stderr, "Bad file creation: %s.\n", filename);
		return -1;
	}

	if(!fprintf(fp, "%2s\n%d %d\n%d\n", pgm->magic, pgm->width, pgm->height, pgm->maximum_value)) {
		fprintf(stderr, "Impossible to write file header.\n");
		fclose(fp);
	}

	int two_bytes = pgm->maximum_value > UINT8_MAX;
	size_t img_elements = pgm->width * pgm->height;
	byte word[2];

	if (two_bytes && ENDIAN_SWAP) {
		for (size_t written_e = 0; written_e < img_elements; ++written_e) {
			//buffer = ENDIAN_ADJUST(((dbyte*)pgm->data)[written_b]);
			word[0] = pgm->data[2 * written_e + 1];
			word[1] = pgm->data[2 * written_e];
			if (!fwrite(word, sizeof(dbyte), 1, fp)) {
				fprintf(stderr, "Impossible to write file data.\n");
				fclose(fp);
				return -1;
			}
			
		}
	} else {
		if (!fwrite(pgm->data, sizeof(byte), img_elements, fp)) {
			fprintf(stderr, "Impossible to write file data.\n");
			fclose(fp);
			return -1;
		}
	}
	
	fclose(fp);

	return 0;
}

void free_pgm(pgm_file* const pgm)
{
	free(pgm->data);
	pgm->data = NULL;
}
