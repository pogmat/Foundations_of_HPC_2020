#include <stdio.h>
#include <stdlib.h>
#include "pgm_bin.h"

#ifndef __GNUC__
#error I need GNU C extensions
#endif

#define CLEANUP(x) __attribute__((__cleanup__(x)))

void cleanup_pgm(pgm_file_t** fpp) {close_pgm(*fpp);}

int main(int argc, char** argv)
{
	if (argc != 3) {
		fprintf(stderr, "USAGE:\n%s <input> <output>\n", argv[0]);
		return EXIT_FAILURE;
	}

	pgm_file_t* CLEANUP(cleanup_pgm) f_in = open_pgm(argv[1], 'r');
	if (!f_in) {
		fprintf(stderr, "Error in input file.\n");
		return EXIT_FAILURE;
	}
       
	pgm_file_t* CLEANUP(cleanup_pgm) f_out = open_pgm(argv[2], 'w');
	if (!f_out) {
		fprintf(stderr, "Error in input file.\n");
		return EXIT_FAILURE;
	}

	if (read_pgm(f_in)) {
		fprintf(stderr, "Error in reading.\n");
		return EXIT_FAILURE;
	}

	// Move the data
	f_out->image = f_in->image;
	f_in->image.data = NULL;

	if (write_pgm(f_out)) {
		fprintf(stderr, "Error in writing.\n");
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
