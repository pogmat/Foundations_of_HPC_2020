/*
 * test.c
 * Simple program to test pgm_bin functionalities.
 *
 * Author: Matteo Poggi <matteo.poggi.fi@gmail.com>
 */ 

#include <stdio.h>
#include "pgm_bin.h"

int main(int argc, char** argv)
{
	if (argc != 3) {
		fprintf(stderr, "Give me two files!\n");
		return 1;
	}

	pgm_file my_pgm;

	if (read_pgm(argv[1], &my_pgm))
		return 1;

	printf("Content of file %s:\n", argv[1]);
	printf("magic value: %s\n", my_pgm.magic);
	printf("height: %d\n", my_pgm.height);
        printf("width:  %d\n", my_pgm.width);
	printf("maximum_value: %d\n", my_pgm.maximum_value);

	if (write_pgm(argv[2], &my_pgm))
		return 1;
	
	return 0;
}       
