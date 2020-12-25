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

void cleanup_pgm_file(pgm_file_t** f) {close_pgm(*f);}
void cleanup_kernel_t(kernel_t** k) {free_kernel(*k);}

int blur_pgm(const pgm_file_t* const input,
	     const kernel_t* const kernel,
	     const pgm_file_t* const output)
{
	int w = (int)input->image.width;
	int h = (int)input->image.height;
	int halo = (int)kernel->s;

	int threads, w_th, h_th;
	#pragma omp parallel firstprivate(w, h, halo)
	{
		#pragma omp master
		{
			// Compute the distribution of w and h threads
			threads = omp_get_num_threads();
			grid_dimension(w, h, threads, &w_th, &h_th);
		}
		
		#pragma omp barrier
		int id = omp_get_thread_num();
		int j = id % w_th;
		int i = id / w_th;
		thimg_t frame;
		// Compute the grid parameters
		get_frame(w, h, w_th, h_th, j, i, halo, &frame);

		#ifdef DEBUG_MODE
		// In debug mote print some info
		#pragma omp critical(show)
		printf("[THREAD %d]  coordinates       : (%d, %d)\n"
		       "            reading dimensions: (%d, %d)\n"
		       "            reading start     : (%d, %d)\n"
		       "            writing dimensions: (%d, %d)\n"
		       "            writing start     : (%d, %d)\n",
		       id, j, i,
		       frame.read.w, frame.read.h, frame.read.j, frame.read.i,
		       frame.writ.w, frame.writ.h, frame.writ.j, frame.writ.i);
		#endif

		// Every thread make a copy of the kernel
		// to minimize memory contention
		kernel_t* CLEANUP(cleanup_kernel_t) private_k = copy_kernel(kernel);
		if (!private_k) {
			fprintf(stderr, "Bad allocation on thread %d.\n", id);
		}

		// Touch first paradigm
		for (int j = frame.writ.j; j < frame.writ.j + frame.writ.w; ++j)
			for (int i = frame.writ.i; i < frame.writ.i + frame.writ.h; ++i)
				touch_zero(&input->image, j + i * w);

		#pragma omp barrier
		#pragma omp master
		if (read_pgm(input)) {
			fprintf(stderr, "Error reading file.\n");
		}
		#pragma omp barrier

		
	}

	return 0;

}

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
	normalize_luminosity(kernel, 1.0);

	pgm_file_t* CLEANUP(cleanup_pgm_file) out_file = open_pgm(argv[3], 'w');
	if (!out_file) {
		fprintf(stderr, "Error in reading output file: %s.\n", argv[3]);
		return EXIT_FAILURE;
	}

	blur_pgm(in_file, kernel, out_file);
	
	return EXIT_SUCCESS;
}
