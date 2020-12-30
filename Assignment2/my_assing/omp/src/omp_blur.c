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

void cleanup_pgm_file(pgm_file_t** f) {close_pgm(f);}
void cleanup_kernel_t(kernel_t** k) {free_kernel(k);}

#ifdef PRIVATE_SUBIMAGE
void cleanup_bytes(byte** b) {free(*b); *b = NULL;}
#endif

int blur_pgm(const pgm_file_t* const input,
	     const kernel_t* const kernel,
	      pgm_file_t* const output)
{
	const int w = (int)input->image.width;
	const int h = (int)input->image.height;
	const int halo = (int)kernel->s;

	int loaded_kernel = 0;
	int loaded_images = 0;
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
			printf("### width  : %8d\n"
			       "### height : %8d\n"
			       "### halo   : %8d\n"
			       "### threads: %8d\n",
			       w, h, halo, threads);
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

		#ifndef NO_TOUCH_FIRST
		// Touch first paradigm
		for (int j = frame.writ.j; j < frame.writ.j + frame.writ.w; ++j)
			for (int i = frame.writ.i; i < frame.writ.i + frame.writ.h; ++i)
				touch_zero(&input->image, j + i * w);
		#else
		#warning NO_TOUCH_FIRST enabled
		#endif
		
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

		#ifndef SHARED_KERNEL
		// Every thread make a copy of the kernel
		// to minimize memory contention
		kernel_t* CLEANUP(cleanup_kernel_t) private_k = copy_kernel(kernel);
		if (!private_k) {
			fprintf(stderr, "[THREAD %d] Bad allocation on thread. (private kernel)\n", id);
		} else {
			#pragma omp atomic update
			++loaded_kernel;
		}
		#else
		#warning SHARED_KERNEL enabled
		#endif

		// ##### END INIT PHASE #####
		T_STOP;
		init_time = DELTA_mSEC;
		
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

			#ifdef SHARED_KERNEL
			loaded_kernel = threads;
			#warning SHARED_KERNEL enabled
			#endif

			#ifndef PRIVATE_SUBIMAGE
			loaded_images = threads;
			#else
			#warning PRIVATE_SUBIMAGE enabled
			#endif
			
			size_t bytes_data = (1 + (input->image.maximum_value > UINT8_MAX))
				* input->image.width * input->image.height;
			if (read_pgm(input)) {
				fprintf(stderr, "Error in reading input file.\n");
			} else if (!(output->image.data = (byte*)malloc(bytes_data))) {
				fprintf(stderr, "Error in allocating output file buffer.\n");
			} else if (loaded_kernel == threads || loaded_images == threads) {
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

		#ifdef PRIVATE_SUBIMAGE
		#warning PRIVATE_SUBIMAGE enabled
		// ##### BEGIN INIT2 PHASE #####
		T_START;

		byte* CLEANUP(cleanup_bytes) private_i = (byte*)malloc(
			frame.read.w * frame.read.h
			* (1 + (input->image.maximum_value > UINT8_MAX))
			* sizeof(byte));

		if (input->image.maximum_value > UINT8_MAX) {
			for (int j = 0; j < frame.read.w; ++j)
				for (int i = 0; i < frame.read.h; ++i)
					((dbyte*)private_i)[j + frame.read.w * i] =
						((dbyte*)(input->image.data))[j + frame.read.j + w * (i + frame.read.i)];
		} else {
			for (int j = 0; j < frame.read.w; ++j)
				for (int i = 0; i < frame.read.h; ++i)
					private_i[j + frame.read.w * i] =
						input->image.data[j + frame.read.j + w * (i + frame.read.i)];
		}
		
		// ##### END INIT2 PHASE #####
		T_STOP;
		init_time += DELTA_mSEC;
		#endif
		init_time2 = init_time * init_time;

		
		// ##### BEGIN COMP PHASE #####
		T_START;
		
		if (proceed) {
			register int s = halo;
			register real new_value;
			#ifndef SHARED_KERNEL
			real* k_p = private_k->kernel;
			#else
			real* k_p = kernel->kernel;
			#warning SHARED KERNEL enabled
			#endif

			#ifdef DEBUG_MODE
			printf("Blurring form kernel %d.\n", id);
			#endif
			
			if (input->image.maximum_value > UINT8_MAX) {
				#ifdef DEBUG_MODE
				#pragma omp master
				printf("16-bit image.\n");
				#endif
				dbyte* n_p = (dbyte*)output->image.data;

				#ifndef PRIVATE_SUBIMAGE
				dbyte* o_p = (dbyte*)input->image.data;
				#else
				#warning PRIVATE_SUBIMAGE enabled
				dbyte* o_p = (dbyte*)private_i;				
				#endif
				
				for (int i = frame.writ.i; i < frame.writ.i + frame.writ.h; ++i)
					for (int j = frame.writ.j; j < frame.writ.j + frame.writ.w; ++j) {
						new_value = 0.0;
						for (int a = max(0, i - s); a < min(h, i + s + 1); ++a) {
							for (int b = max(0, j - s); b < min(w, j + s +1) - 1; b += 2) {
								#ifndef PRIVATE_SUBIMAGE
								new_value +=
									(o_p[b + w * a] * k_p[b - j + s + (2 * s +1) * (a - i + s)] +
									 o_p[b + 1 + w * a] * k_p[b + 1 - j + s + (2 * s +1) * (a - i + s)]);
								#else
								#warning PRIVATE_SUBIMAGE enabled
								new_value +=
									(o_p[b - frame.read.j + frame.read.w * (a - frame.read.i)] * k_p[b - j + s + (2 * s +1) * (a - i + s)] +
									 o_p[b - frame.read.j + 1 + frame.read.w * (a - frame.read.i)] * k_p[b + 1 - j + s + (2 * s +1) * (a - i + s)]);
								#endif
							}
							if ((min(w, j + s +1) - max(0, j - s)) % 2) {
								#ifndef PRIVATE_SUBIMAGE
								new_value +=
									o_p[min(w, j + s +1) - 1 + w * a] * k_p[min(w, j + s +1) - 1 - j + s + (2 * s +1) * (a - i + s)];
								#else
								#warning PRIVATE_SUBIMAGE enabled
								new_value +=
									o_p[min(w, j + s +1) - 1 - frame.read.j+ frame.read.w * (a - frame.read.i)]
									* k_p[min(w, j + s +1) - 1 - j + s + (2 * s +1) * (a - i + s)];							       
								#endif
							}
						}
						n_p[j + w * i] = (dbyte)min(UINT16_MAX, (uint64_t)(new_value + 0.5));
						//printf("### %d\n", n_p[j + w * i]);
					}

			} else {
				#ifdef DEBUG_MODE
				#pragma omp master
				printf("8-bit image.\n");
				#endif
				byte* n_p = output->image.data;

				#ifndef PRIVATE_SUBIMAGE
				byte* o_p = input->image.data;
				#else
				#warning PRIVATE SUBIMAGE enabled
				byte* o_p = private_i;
				#endif

				for (int i = frame.writ.i; i < frame.writ.i + frame.writ.h; ++i)
					for (int j = frame.writ.j; j < frame.writ.j + frame.writ.w; ++j) {
						new_value = 0.0;
						for (int a = max(0, i - s); a < min(h, i + s + 1); ++a) {
							for (int b = max(0, j - s); b < min(w, j + s +1) - 1; b += 2) {
								#ifndef PRIVATE_SUBIMAGE
								new_value +=
									(o_p[b + w * a] * k_p[b - j + s + (2 * s +1) * (a - i + s)] +
									 o_p[b + 1 + w * a] * k_p[b + 1 - j + s + (2 * s +1) * (a - i + s)]);
								#else
								#warning PRIVATE SUBIMAGE enabled
								new_value +=
									(o_p[b - frame.read.j + frame.read.w * (a - frame.read.i)] * k_p[b - j + s + (2 * s +1) * (a - i + s)] +
									 o_p[b - frame.read.j + 1 + frame.read.w * (a - frame.read.i)] * k_p[b + 1 - j + s + (2 * s +1) * (a - i + s)]);
								#endif
							}
							if ((min(w, j + s +1) - max(0, j - s)) % 2) {
								#ifndef PRIVATE_SUBIMAGE
								new_value +=
									o_p[min(w, j + s +1) - 1 + w * a] * k_p[min(w, j + s +1) - 1 - j + s + (2 * s +1) * (a - i + s)];
								#else
								#warning PRIVATE SUBIMAGE enabled
								new_value +=
									o_p[min(w, j + s +1) - 1 - frame.read.j+ frame.read.w * (a - frame.read.i)]
									* k_p[min(w, j + s +1) - 1 - j + s + (2 * s +1) * (a - i + s)];							       
								#endif								
							}
						}
						n_p[j + w * i] = (byte)min(UINT8_MAX, (uint64_t)(new_value + 0.5));
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

	printf("init_time: %14.6lf +- %11.6lf ms\n", init_mean, init_stddev);
	printf("load_time: %14.6lf                ms\n", load_time);
	printf("comp_time: %14.6lf +- %11.6lf ms\n", comp_mean, comp_stddev);
	printf("flsh_time: %14.6lf                ms\n", flsh_time);

	if (!proceed)
		return 1;

	return 0;

}

int main(int argc, char** argv)
{
	double tot_time = 0;
	struct timespec TIMESTAMP_START, TIMESTAMP_STOP;

	T_START;
	
	if (argc < 4) {
	usage:
		fprintf(stderr, "USAGE:\n%s [kernel-type] [kernel-size] {kernel-param} [input-file] {output-file}\n", argv[0]);
		return EXIT_FAILURE;
	}

	char default_out[] = "out.pgm";
	
	int kernel_type = atoi(argv[1]);
	unsigned int kernel_size = (unsigned int)atoi(argv[2]);
        char* output_file = default_out;
	char* input_file;
	double kernel_param;
	
	switch (kernel_type) {
	case 0:
	case 2:
		input_file = argv[3];
		if (argc == 5)
			output_file = argv[4];
		if (argc > 5)
			goto usage;
		break;
	case 1:
		if (argc < 5) {
			fprintf(stderr, "Weighted kernel needs one more parameter.\n");
			return EXIT_FAILURE;
		}
		kernel_param = atof(argv[3]);
		input_file = argv[4];
		if (argc == 6)
			output_file = argv[5];
		if (argc > 6)
			goto usage;
		break;
	default:
		fprintf(stderr, "Unkown kernel:\n1. Mean kernel\n2. Weighted kernel\n3. Gaussian kernel\n");
		return EXIT_FAILURE;
	}

	if (!(kernel_size % 2)) {
		fprintf(stderr, "Kernel size must be odd.\n");
		return EXIT_FAILURE;
	}

	int s = kernel_size / 2;

	pgm_file_t* CLEANUP(cleanup_pgm_file) in_file = open_pgm(input_file, 'r');
	if (!in_file) {
		fprintf(stderr, "Error in reading input file: %s.\n", argv[1]);
		return EXIT_FAILURE;
	}
	
	kernel_t* CLEANUP(cleanup_kernel_t) kernel = NULL;

	switch (kernel_type) {
	case 0:
		kernel = mean_kernel((unsigned int)s);
		break;
	case 1:
		kernel = weight_kernel((unsigned int)s, kernel_param);
		break;
	case 2:
		kernel = gaussian_kernel((unsigned int)s);
	}

	if (!kernel) {
		fprintf(stderr, "Error in reading kernel file: %s.\n", argv[2]);
		return EXIT_FAILURE;
	}
	
	pgm_file_t* CLEANUP(cleanup_pgm_file) out_file = open_pgm(output_file, 'w');
	if (!out_file) {
		fprintf(stderr, "Error in reading output file: %s.\n", argv[3]);
		return EXIT_FAILURE;
	}

	if (blur_pgm(in_file, kernel, out_file)) {
		fprintf(stderr, "Error in blurring image.\n");
		return EXIT_FAILURE;
	}
	
	
	/*
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
	*/

	T_STOP;
	tot_time = DELTA_mSEC;
	printf("totl_time: %14.6lf                ms\n", tot_time);
	
	return EXIT_SUCCESS;
}
