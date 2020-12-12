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
	size_t bytes_data;
	size_t read_b;
	
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
	
	bytes_data = ((pgm->maximum_value > 255) ? 2 : 1) * pgm->width * pgm->height;

	pgm->data = (byte*)malloc(bytes_data * sizeof(byte));
	if (!pgm->data) {
		fprintf(stderr, "Bad allocation.\n");
		fclose(fp);
		return -1;
	}
	
	read_b = fread(pgm->data, sizeof(byte), bytes_data, fp);
	if (read_b != bytes_data) {
		fprintf(stderr, "Invalid format: the file is corrupted.\n");
		fclose(fp);
		return -1;
	}

	fclose(fp);
	
	return 0;
}

int write_pgm(const char* const filename, const pgm_file* const pgm)
{
	FILE* fp = fopen(filename, "wb");
	unsigned long bytes_data;
	size_t written_b;

	if (fp == NULL) {
		fprintf(stderr, "Bad file creation: %s.\n", filename);
		return -1;
	}

	if(!fprintf(fp, "%2s\n%d %d\n%d\n", pgm->magic, pgm->width, pgm->height, pgm->maximum_value)) {
		fprintf(stderr, "Impossible to write file header.\n");
		fclose(fp);
	}

	bytes_data = ((pgm->maximum_value > 255) ? 2 : 1) * pgm->width * pgm->height;

	written_b = fwrite(pgm->data, sizeof(byte), bytes_data, fp);
	if (written_b != bytes_data) {
		fprintf(stderr, "Impossible to write file data.\n");
		fclose(fp);
		return -1;
	}
		
	fclose(fp);

	return 0;
}

void free_pgm(pgm_file* const pgm)
{
	free(pgm->data);
	pgm->data = NULL;
}
