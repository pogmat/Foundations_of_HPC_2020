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
	size_t img_size = (two_bytes ? 2 : 1) * pgm->width * pgm->height;
	byte word[2];

	pgm->data = (byte*)malloc(img_size);
	if (!pgm->data) {
		fprintf(stderr, "Bad allocation.\n");
		fclose(fp);
		return -1;
	}
	
	for (size_t read_b = 0; read_b < (img_size / 2) * 2; read_b += 2) {
		if (!(fread(word, 2, 1, fp))) {
			fprintf(stderr, "Invalid format: the file is corrupted\n");
			fclose(fp);
			return -1;
		}

		pgm->data[read_b] = two_bytes ? word[1] : word[0];
		pgm->data[read_b + 1] = two_bytes ? word[0] : word[1];
	}

	if (img_size % 2) {
		if (!fread(pgm->data + img_size - 1, 1, 1, fp)) {
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
	size_t img_size = (two_bytes ? 2 : 1) * pgm->width * pgm->height;
	byte word[2];

	for (size_t written_b = 0; written_b < (img_size / 2) * 2; written_b += 2) {
		word[0] = two_bytes ? pgm->data[written_b + 1] : pgm->data[written_b];
		word[1] = two_bytes ? pgm->data[written_b] : pgm->data[written_b + 1];

		if (fwrite(word, 1, 2, fp) != 2) {
			fprintf(stderr, "Impossible to write file data.\n");
			fclose(fp);
			return -1;
		}
	}

	if (img_size % 2) {
		if (!fwrite(pgm->data + img_size -1, 1, 1, fp)) {
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
