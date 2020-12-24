#include <stdio.h>
#include <stdlib.h>
#include "pgm_bin.h"
#include "kernel.h"
#include "grid.h"
#include <omp.h>

#define _GNU_SOURCE

#ifndef __GNUC__
#error This program require GNU C extensions
#endif
#define CLEANUP(x) __attribute__((__cleanup__(x)))

#ifndef _OPENMP
#warning You are compiling without openmp support!
#endif

#ifdef DEBUG_MODE
void test_grid(int w, int h, int k)
{
	int threads;
	#pragma omp parallel
	#pragma omp master
	{
		threads = omp_get_num_threads();
	}

	printf("Operating with %d threads\n", threads);
	printf("Image dimension: %d x %d; halo of %d\n", w, h, k);

	int w_th, h_th;
	grid_dimension(w, h, threads, &w_th, &h_th);
	printf("Distribution of threads: %d x %d\n", w_th, h_th);

	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int j = id % w_th;
		int i = id / w_th;
		thimg_t frame;
		get_frame(w, h, w_th, h_th, j, i, k, &frame);

		#pragma omp critical(show)
		{
			printf("[THREAD %d]  coordinates       : (%d, %d)\n"
			       "            reading dimensions: (%d, %d)\n"
			       "            reading start     : (%d, %d)\n"
			       "            writing dimensions: (%d, %d)\n"
			       "            writing start     : (%d, %d)\n",
			       id, j, i,
			       frame.read.w, frame.read.h, frame.read.j, frame.read.i,
			       frame.writ.w, frame.writ.h, frame.writ.j, frame.writ.i);
		}
	}
}
#endif

void cleanup_pgm_file(pgm_file_t** f) {close_pgm(*f);}
void cleanup_kernel_t(kernel_t** k) {free_kernel(*k);}

int main(int argc, char** argv)
{
	if (argc != 4) {
		fprintf(stderr, "USAGE:\n%s <input> <kernel> <output>\n", argv[0]);
		return EXIT_FAILURE;
	}

	pgm_file_t* CLEANUP(cleanup_pgm_file) in_file = open_pgm(argv[1], 'r');
	if (!in_file) {
		fprintf(stderr, "Error in reading input file: %s.\n", argv[1]);
		return EXIT_FAILURE;
	}

	kernel_t* CLEANUP(cleanup_kernel_t) kernel = init_kernel_from_file(argv[2]);
	if (!kernel) {
		fprintf(stderr, "Error in reading kernel file: %s.\n", argv[2]);
		return EXIT_FAILURE;
	}

	pgm_file_t* CLEANUP(cleanup_pgm_file) out_file = open_pgm(argv[3], 'w');
	if (!out_file) {
		fprintf(stderr, "Error in reading output file: %s.\n", argv[3]);
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
