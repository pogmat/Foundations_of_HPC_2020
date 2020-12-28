#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "pgm_bin.h"
#include "kernel.h"
#include "grid.h"

#define _GNU_SOURCE

#ifndef __GNUC__
#error This program require GNU C extensions
#endif
#define CLEANUP(x) __attribute__((__cleanup__(x)))

#define TIMESTAMP_START t_start
#define TIMESTAMP_STOP t_stop
#define T_START clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(TIMESTAMP_START))
#define T_STOP clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(TIMESTAMP_STOP))
#define DELTA_mSEC 1e3 * (double)(TIMESTAMP_STOP.tv_sec - TIMESTAMP_START.tv_sec) + 1e-6 * (double)(TIMESTAMP_STOP.tv_nsec - TIMESTAMP_START.tv_nsec) 

#ifndef _OPENMP
#warning You are compiling without openmp support!
#endif

void cleanup_pgm_file(pgm_file_t** f) {close_pgm(*f);}
void cleanup_kernel_t(kernel_t** k) {free_kernel(*k);}

int blur_pgm(const pgm_file_t* const input,
	     const kernel_t* const kernel,
	      pgm_file_t* const output)
{
	const int w = (int)input->image.width;
	const int h = (int)input->image.height;
	const int halo = (int)kernel->s;
	int loaded_kernel = 0;
	int proceed = 0;
       	int threads, w_th, h_th;
	double init_time = 0,
		init_time2 = 0,
		load_time = 0,
		comp_time = 0,
		comp_time2 = 0,
		flsh_time = 0;
	
	#pragma omp parallel firstprivate(w, h, halo) reduction(+:init_time, init_time2, comp_time, comp_time2)
	{
		struct timespec TIMESTAMP_START, TIMESTAMP_STOP;
		#pragma omp master
		{
			// Compute the distribution of w and h threads
			threads = omp_get_num_threads();
			grid_dimension(w, h, threads, &w_th, &h_th);
		}
		#pragma omp barrier

		// ##### BEGIN INIT PHASE #####
		T_START;
		
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
		
		// Touch first paradigm
		for (int j = frame.writ.j; j < frame.writ.j + frame.writ.w; ++j)
			for (int i = frame.writ.i; i < frame.writ.i + frame.writ.h; ++i)
				touch_zero(&input->image, j + i * w);

		#ifdef DEBUG_MODE
		#pragma omp master
		{
			printf("First touch done.\n");
			printf("KERNEL:\n");
			for (int j = 0; j < (2*halo+1); ++j) {
				for (int i = 0; i < (2*halo+1); ++i)
					printf("%1.5lf  ", kernel->kernel[j + i * (2*halo+1)]);
				printf("\n");
			}
		}
		#endif
			    
		// Every thread make a copy of the kernel
		// to minimize memory contention
		kernel_t* CLEANUP(cleanup_kernel_t) private_k = copy_kernel(kernel);
		if (!private_k) {
			fprintf(stderr, "[THREAD %d] Bad allocation on thread.\n", id);
		} else {
			#pragma omp atomic update
			++loaded_kernel;
		}

		// ##### END INIT PHASE #####
		T_STOP;
		init_time = DELTA_mSEC;
		init_time2 = init_time * init_time;
		
		// Data are written from file to input buffer
		// Memory is allocated for output buffer
		#pragma omp barrier
		#pragma omp master
	        {

			// ##### BEGIN LOAD PHASE #####
			T_START;
			
			#ifdef DEBUG_MODE
			printf("Private kernels: %d/%d loaded.\n", loaded_kernel, threads);
			#endif
			
			size_t bytes_data = (1 + (input->image.maximum_value > UINT8_MAX))
				* input->image.width * input->image.height;
			if (read_pgm(input)) {
				fprintf(stderr, "Error in reading input file.\n");
			} else if (!(output->image.data = (byte*)malloc(bytes_data))) {
				fprintf(stderr, "Error in allocating output file buffer.\n");
			} else if (loaded_kernel == threads) {
				strncpy(output->image.magic, input->image.magic, 3);
				output->image.width = input->image.width;
				output->image.height = input->image.height;
				output->image.maximum_value = input->image.maximum_value;
				proceed = 1;
			}

			#ifdef DEBUG_MODE
			printf("proceed = %d.\n", proceed);
			#endif

			// ##### END LOAD PHASE #####
			T_STOP;
			load_time = DELTA_mSEC;
		}		
		#pragma omp barrier

		// ##### BEGIN COMP PHASE #####
		T_START;
		
		if (proceed) {
			register int s = halo;
			register real new_value;
			real* k_p = private_k->kernel;

			#ifdef DEBUG_MODE
			printf("Blurring form kernel %d.\n", id);
			#endif
			
			if (input->image.maximum_value > UINT8_MAX) {
				#ifdef DEBUG_MODE
				#pragma omp master
				printf("16-bit image.\n");
				#endif
				dbyte* n_p = (dbyte*)output->image.data;
				dbyte* o_p = (dbyte*)input->image.data;

				for (int i = frame.writ.i; i < frame.writ.i + frame.writ.h; ++i)
					for (int j = frame.writ.j; j < frame.writ.j + frame.writ.w; ++j) {
						new_value = 0.0;
						for (int a = max(0, i - s); a < min(h, i + s + 1); ++a) {
							for (int b = max(0, j - s); b < min(w, j + s +1) - 1; b += 2) {
								new_value +=
									(o_p[b + w * a] * k_p[b - j + s + (2 * s +1) * (a - i + s)] +
									 o_p[b + 1 + w * a] * k_p[b + 1 - j + s + (2 * s +1) * (a - i + s)]);
							}
							if ((min(w, j + s +1) - max(0, j - s)) % 2)
								new_value +=
									o_p[min(w, j + s +1) - 1 + w * a] * k_p[min(w, j + s +1) - 1 - j + s + (2 * s +1) * (a - i + s)];
						}
						n_p[j + w * i] = (dbyte)min(UINT16_MAX, (uint64_t)(new_value + 0.5));
					}
			} else {
				#ifdef DEBUG_MODE
				#pragma omp master
				printf("8-bit image.\n");
				#endif
				byte* n_p = output->image.data;
				byte* o_p = input->image.data;
				
				for (int i = frame.writ.i; i < frame.writ.i + frame.writ.h; ++i)
					for (int j = frame.writ.j; j < frame.writ.j + frame.writ.w; ++j) {
						new_value = 0.0;
						for (int a = max(0, i - s); a < min(h, i + s + 1); ++a) {
							for (int b = max(0, j - s); b < min(w, j + s +1) - 1; b += 2) {
								new_value +=
									(o_p[b + w * a] * k_p[b - j + s + (2 * s +1) * (a - i + s)] +
									 o_p[b + 1 + w * a] * k_p[b + 1 - j + s + (2 * s +1) * (a - i + s)]);						
							}
							if ((min(w, j + s +1) - max(0, j - s)) % 2)
								new_value +=
									o_p[min(w, j + s +1) - 1 + w * a] * k_p[min(w, j + s +1) - 1 - j + s + (2 * s +1) * (a - i + s)];
						}
						n_p[j + w * i] = (dbyte)min(UINT8_MAX, (uint64_t)(new_value + 0.5));
					}								
			}

			// ##### END COMP PHASE #####
			T_STOP;
			comp_time = DELTA_mSEC;
			comp_time2 = comp_time * comp_time;			
			
			#pragma omp barrier
			#pragma omp master
			{

				// ##### BEGIN FLSH PHASE #####
				T_START;
				
				#ifdef DEBUG_MODE
				printf("Blurred.\n");
				#endif
				
				if (write_pgm(output)) {
					fprintf(stderr, "Error in file writing.\n");
					proceed = 0;
				}
				
				#ifdef DEBUG_MODE
				printf("Image written.\n");
				#endif

				// ##### END FLSH PHASE #####
				T_STOP;
				flsh_time = DELTA_mSEC;
			}
		}
	}

	double init_mean = init_time / threads;
	double init_stddev = sqrt(init_time2 / threads - init_mean * init_mean);
	double comp_mean = comp_time / threads;
	double comp_stddev = sqrt(comp_time2 / threads - comp_mean * comp_mean);

	printf("init_time: %11.6lf +- %11.6lf ms\n", init_mean, init_stddev);
	printf("load_time: %11.6lf                ms\n", load_time);
	printf("comp_time: %11.6lf +- %11.6lf ms\n", comp_mean, comp_stddev);
	printf("load_time: %11.6lf                ms\n", flsh_time);

	if (!proceed)
		return 1;

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

	if (blur_pgm(in_file, kernel, out_file)) {
		fprintf(stderr, "Error in blurring image.\n");
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
