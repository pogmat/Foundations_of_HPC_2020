#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#ifndef __GNUC__
#error This program relies on GNU C extension
#endif

#define DEBUG_VAR(x) printf("DEBUG "#x": %d\n", x)

#define TRIALS 15
#define MY_CLOCK CLOCK_REALTIME
#define CLEANUP(x) __attribute__((__cleanup__(x)))


typedef struct {
	int w;			// image width
	int h;			// imshr height
} image_t;

typedef struct {
	const image_t* image_p;
	int t_w;		// columns of the grid
	int t_h;		// rows of the grid
} thread_image_t;

typedef struct {
	const thread_image_t* t_image_p;
	int rw;			// width of the frame (reading)
	int rh;			// height of the frame (reading)
	int ir_start;		// starting i (reading)
	int jr_start;		// starting j (reading)
	int ww;			// width of the frame (writing)
	int wh;			// height of the frame (writing)
	int iw_start;		// starting i (writing)
	int jw_start;		// starting j (writing)
} framed_t;

static inline int min(const int a, const int b) {return((a < b) ? a : b);}
static inline int max(const int a, const int b) {return((a > b) ? a : b);}

static void cleanup_uint16(uint16_t** pp)
{
	if (*pp)
		free(*pp);
}

static void cleanup_double(double** pp)
{
	if (*pp)
		free(*pp);
}

void grid_dimension(const image_t* const img, const int threads, thread_image_t* const t_img)
{

	int big_d, * big_p, * small_p;

	if (img->w > img->h) {
		big_d = img->w;
		big_p = &t_img->t_w;
		small_p = &t_img->t_h;
	} else {
		big_d = img->h;
		big_p = &t_img->t_h;
		small_p = &t_img->t_w;
	}
	
	int i;
	i = (((float)big_d / (sqrtf(img->h * img->w / threads))) + 0.5);
	
	while (1) {
		if (threads % i)
			++i;
		else {
			*big_p = i;
			*small_p = threads / i;
			break;
		}	
	}

	t_img->image_p = img;
}

void get_frame(const thread_image_t* const t_img, int halo, int i, int j, framed_t* const frame)
{
	frame->t_image_p = t_img;

	int w = t_img->image_p->w;
	int h = t_img->image_p->h;
	int t_w = t_img->t_w;
	int t_h = t_img->t_h;

	int chunk_w = w / t_w;
	int chunk_h = h / t_h;
	
	int thed_w = t_w - w + chunk_w * t_w - 1;
	int thed_h = t_h - h + chunk_h * t_h - 1;

	int ww = chunk_w + (j > thed_w);
	int wh = chunk_h + (i > thed_h);
	frame->ww = ww;
	frame->wh = wh;
	
	int tl_j = chunk_w * j + max(0, j - thed_w - 1);
	int tl_i = chunk_h * i + max(0, i - thed_h - 1);
	int br_j = tl_j + ww - 1;
	int br_i = tl_i + wh - 1;

	frame->iw_start = tl_i;
	frame->jw_start = tl_j;
	
	int offset_l = min(halo, tl_j);
	int offset_r = min(halo, w - br_j - 1);
	int offset_t = min(halo, tl_i);
	int offset_b = min(halo, h - br_i - 1);
	
	ww += offset_l + offset_r;
	wh += offset_t + offset_b;

	frame->rw = ww;
	frame->rh = wh;
	frame->ir_start = tl_i - offset_t;
	frame->jr_start = tl_j - offset_l;
	
}

int touch_by_one(int w, int h, int k)
{
	image_t img = {.w = w, .h = h};
	long int kernel_dim = (2 * k + 1) * (2 * k + 1);

	// Calloc here. This is touch by one.
	uint16_t* input CLEANUP(cleanup_uint16) = (uint16_t*)calloc(w * h, sizeof(uint16_t));
	if (!input) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}
	// Calloc here. This is touch by one.
	uint16_t* output CLEANUP(cleanup_uint16) = (uint16_t*)calloc(w * h, sizeof(uint16_t));
	if (!output) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}
	// Calloc here. This is touch by one.
	double* kernel CLEANUP(cleanup_double) = (double*)calloc(kernel_dim, sizeof(double));
	if (!kernel) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}	
	
	long pseudo_seed = (long)(strlen(getenv("PATH")) * getpid())^(long)(&img);
	struct drand48_data rnd;
	double pseudo_random;
	srand48_r(pseudo_seed, &rnd);

	for (int i = 0; i < w * h; ++i)
		input[i] = (drand48_r(&rnd, &pseudo_random),(uint16_t)(pseudo_random * UINT16_MAX + 0.5));

	for (int i = 0; i < kernel_dim; ++i)
		kernel[i] = (drand48_r(&rnd, &pseudo_random), pseudo_random / kernel_dim);
	
	//printf("Data initialized by master.\n");

	int threads;
	thread_image_t t_img;

	
	#pragma omp parallel
	{
		#pragma omp master
		{
			threads = omp_get_num_threads();
			//printf("threads number: %d\n", threads);
			//printf("allocated memory: %f Mb.\n", ((float)(w*h*sizeof(int16_t))) / ((float)(2 << 19)));
			grid_dimension(&img, threads, &t_img);
			//printf("w_threads: %d, h_threads: %d.\n\n", t_img.t_w, t_img.t_h);	
		}
		
		#pragma omp barrier
		
		
		int myid = omp_get_thread_num();
		int i = myid / t_img.t_w;
		int j = myid % t_img.t_w;
		framed_t frame;
		//get_frame(&t_img, k, i, j, &frame);
		get_frame(&t_img, k, i, j, &frame);
		
		//#pragma omp critical(show)
		//{
			//printf("[THREAD %d]  coordinates       : (%d, %d)\n"
			//       "            reading dimensions: (%d, %d)\n"
			//       "            reading start     : (%d, %d)\n"
			//       "            writing dimensions: (%d, %d)\n"
			//       "            writing start     : (%d, %d)\n",
			//       myid, j, i, frame.rw, frame.rh, frame.jr_start, frame.ir_start,
			//       frame.ww, frame.wh, frame.jw_start, frame.iw_start);
			//}
		
		double new_value;

		for (int i = frame.iw_start; i < frame.iw_start + frame.wh; ++i)
			for (int j = frame.jw_start; j < frame.jw_start + frame.ww; ++j) {
				new_value = 0.0;
				for (int a = max(0, i - k); a < min(frame.wh, i + k + 1); ++a)
					for (int b = max(0, j - k); b < min(frame.ww, j + k + 1); ++b)
						new_value += input[b + frame.ww * a]
							* kernel[b - j + k + (2 * k + 1) * (a - i + k)];
				output[j + frame.ww * i] = (uint16_t)min(UINT16_MAX, (int)(new_value + 0.5));
			}
	}
	
	// Prevent non-computation of output.
	double index;
	drand48_r(&rnd, &index);
	volatile uint16_t sink = output[w * h * (int)index];
	(void) sink;
	
	return 0;
}

int touch_by_all_common_kernel(int w, int h, int k)
{
	image_t img = {.w = w, .h = h};
	long int kernel_dim = (2 * k + 1) * (2 * k + 1);

	// Malloc here. Memory is only allocated and not initialized.
	uint16_t* input CLEANUP(cleanup_uint16) = (uint16_t*)malloc(w * h * sizeof(uint16_t));
	if (!input) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}
	// Malloc here. Memory is only allocated and not initialized.
	uint16_t* output CLEANUP(cleanup_uint16) = (uint16_t*)malloc(w * h * sizeof(uint16_t));
	if (!output) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}
	// Calloc here. This is touch by one.
	double* kernel CLEANUP(cleanup_double) = (double*)calloc(kernel_dim, sizeof(double));
	if (!kernel) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}	
	
	long pseudo_seed = (long)(strlen(getenv("PATH")) * getpid())^(long)(&img);
	struct drand48_data rnd;
	double pseudo_random;
	srand48_r(pseudo_seed, &rnd);	

	int threads;
	thread_image_t t_img;
	
	#pragma omp parallel
	{
		#pragma omp master
		{
			threads = omp_get_num_threads();
			grid_dimension(&img, threads, &t_img);
		}

		#pragma omp barrier
		
		int myid = omp_get_thread_num();
		int i = myid / t_img.t_w;
		int j = myid % t_img.t_w;
		framed_t frame;
		get_frame(&t_img, k, i, j, &frame);

		// This is touch first.
		for (int i = frame.iw_start; i < frame.iw_start + frame.wh; ++i)
			for (int j = frame.jw_start; j < frame.jw_start + frame.ww; ++j) {
				input[j + frame.ww * i] = 0;
			}

		//Perhaps it would be better to have sigle here.
		//However one should change the initialization above.
		#pragma omp master
		{
			for (int i = 0; i < w * h; ++i)
				input[i] = (drand48_r(&rnd, &pseudo_random),(uint16_t)(pseudo_random * UINT16_MAX + 0.5));

			for (int i = 0; i < kernel_dim; ++i)
				kernel[i] = (drand48_r(&rnd, &pseudo_random), pseudo_random / kernel_dim);
		}

		#pragma omp barrier
				
		double new_value;
		
		for (int i = frame.iw_start; i < frame.iw_start + frame.wh; ++i)
			for (int j = frame.jw_start; j < frame.jw_start + frame.ww; ++j) {
				new_value = 0.0;
				for (int a = max(0, i - k); a < min(frame.wh, i + k + 1); ++a)
					for (int b = max(0, j - k); b < min(frame.ww, j + k + 1); ++b)
						new_value += input[b + frame.ww * a]
							* kernel[b - j + k + (2 * k + 1) * (a - i + k)];
				output[j + frame.ww * i] = (uint16_t)min(UINT16_MAX, (int)(new_value + 0.5));
			}
	}
	
	// Prevent non-computation of output.
	double index;
	drand48_r(&rnd, &index);
	volatile uint16_t sink = output[w * h * (int)index];
	(void) sink;
	
	return 0;
		
}

int touch_by_all_copy_kernel(int w, int h, int k)
{
	image_t img = {.w = w, .h = h};
	long int kernel_dim = (2 * k + 1) * (2 * k + 1);

	// Malloc here. Memory is only allocated and not initialized.
	uint16_t* input CLEANUP(cleanup_uint16) = (uint16_t*)malloc(w * h * sizeof(uint16_t));
	if (!input) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}
	// Malloc here. Memory is only allocated and not initialized.
	uint16_t* output CLEANUP(cleanup_uint16) = (uint16_t*)malloc(w * h * sizeof(uint16_t));
	if (!output) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}
	// Calloc here. This is touch by one.
	double* kernel CLEANUP(cleanup_double) = (double*)calloc(kernel_dim, sizeof(double));
	if (!kernel) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}	
	
	long pseudo_seed = (long)(strlen(getenv("PATH")) * getpid())^(long)(&img);
	struct drand48_data rnd;
	double pseudo_random;
	srand48_r(pseudo_seed, &rnd);	

	int threads;
	thread_image_t t_img;
	
	#pragma omp parallel
	{
		#pragma omp master
		{
			threads = omp_get_num_threads();
			grid_dimension(&img, threads, &t_img);
		}

		#pragma omp barrier
		
		int myid = omp_get_thread_num();
		int i = myid / t_img.t_w;
		int j = myid % t_img.t_w;
		framed_t frame;
		get_frame(&t_img, k, i, j, &frame);

		// This is touch first.
		for (int i = frame.iw_start; i < frame.iw_start + frame.wh; ++i)
			for (int j = frame.jw_start; j < frame.jw_start + frame.ww; ++j) {
				input[j + frame.ww * i] = 0;
			}

		// Perhaps it would be better to have sigle here.
		// However one should change the initialization above.
		#pragma omp master
		{
			for (int i = 0; i < w * h; ++i)
				input[i] = (drand48_r(&rnd, &pseudo_random),(uint16_t)(pseudo_random * UINT16_MAX + 0.5));

			for (int i = 0; i < kernel_dim; ++i)
				kernel[i] = (drand48_r(&rnd, &pseudo_random), pseudo_random / kernel_dim);
		}

		#pragma omp barrier

		// This make a local copy of the kernel
		double* kernel_private CLEANUP(cleanup_double) = (double*)malloc(kernel_dim * sizeof(double));
		for (int i = 0; i < kernel_dim; ++i)
			kernel_private[i] = kernel[i];
			
		double new_value;
		
		for (int i = frame.iw_start; i < frame.iw_start + frame.wh; ++i)
			for (int j = frame.jw_start; j < frame.jw_start + frame.ww; ++j) {
				new_value = 0.0;
				for (int a = max(0, i - k); a < min(frame.wh, i + k + 1); ++a)
					for (int b = max(0, j - k); b < min(frame.ww, j + k + 1); ++b)
						new_value += input[b + frame.ww * a]
							* kernel_private[b - j + k + (2 * k + 1) * (a - i + k)];
				output[j + frame.ww * i] = (uint16_t)min(UINT16_MAX, (int)(new_value + 0.5));
			}
	}
	
	// Prevent non-computation of output.
	double index;
	drand48_r(&rnd, &index);
	volatile uint16_t sink = output[w * h * (int)index];
	(void) sink;
	
	return 0;
		
}

int main(int argc, char** argv)
{      	
	if (argc != 4) {
		fprintf(stderr,
			"USAGE:\n%s <width> <height> <kernel_size>\n", argv[0]);

		return 1;
	}

	int w = atoi(argv[1]);
	int h = atoi(argv[2]);
	int k = atoi(argv[3]);

	struct timespec t_start, t_stop;
	double delta_t, sum_t, sum_t2, mean, stdev;

	printf("Touch by one\n");
	sum_t = 0;
	sum_t2 = 0;
	for (int i = 0; i < TRIALS; ++i) {

		clock_gettime(MY_CLOCK, &t_start);
		
		if (touch_by_one(w, h, k)) {
			fprintf(stderr, "Error in touch_by_one\n");
			return 1;
		}

		clock_gettime(MY_CLOCK, &t_stop);
		delta_t = 1e3 * (double)(t_stop.tv_sec - t_start.tv_sec) + 1e-6 * (double)(t_stop.tv_nsec - t_start.tv_nsec);
		//printf("t = %lf\n", delta_t);
		sum_t += delta_t;
		sum_t2 += delta_t * delta_t;
	}

	mean = sum_t / TRIALS;
	stdev = sqrt(sum_t2 / TRIALS - mean * mean);
	printf("------------------------------------");
	printf("t = %lf +- %lf\n", mean, stdev);

	printf("Touch by all (common kernel)\n");
	sum_t = 0;
	sum_t2 = 0;
	for (int i = 0; i < TRIALS; ++i) {

		clock_gettime(MY_CLOCK, &t_start);
		
		if (touch_by_all_common_kernel(w, h, k)) {
			fprintf(stderr, "Error in touch_by_all_common_kernel\n");
			return 1;
		}

		clock_gettime(MY_CLOCK, &t_stop);
		delta_t = 1e3 * (double)(t_stop.tv_sec - t_start.tv_sec) + 1e-6 * (double)(t_stop.tv_nsec - t_start.tv_nsec);
		//printf("t = %lf\n", delta_t);
		sum_t += delta_t;
		sum_t2 += delta_t * delta_t;
	}

	mean = sum_t / TRIALS;
	stdev = sqrt(sum_t2 / TRIALS - mean * mean);
	printf("------------------------------------");
	printf("t = %lf +- %lf\n", mean, stdev);	

	printf("Touch by all (copy kernel)\n");
	sum_t = 0;
	sum_t2 = 0;
	for (int i = 0; i < TRIALS; ++i) {

		clock_gettime(MY_CLOCK, &t_start);
		
		if (touch_by_all_copy_kernel(w, h, k)) {
			fprintf(stderr, "Error in touch_by_all_common_kernel\n");
			return 1;
		}

		clock_gettime(MY_CLOCK, &t_stop);
		delta_t = 1e3 * (t_stop.tv_sec - t_start.tv_sec) + 1e-6 * (t_stop.tv_nsec - t_start.tv_nsec);
		//printf("t = %lf\n", delta_t);
		sum_t += delta_t;
		sum_t2 += delta_t * delta_t;
	}

	mean = sum_t / TRIALS;
	stdev = sqrt(sum_t2 / TRIALS - mean * mean);
	printf("------------------------------------");
	printf("t = %lf +- %lf\n", mean, stdev);
	
	return 0;
}
