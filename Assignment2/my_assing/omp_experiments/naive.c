#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define DEBUG_VAR(x) printf("DEBUG "#x": %d\n", x)

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
	int sw;			// width of the frame
	int sh;			// height of the frame
	int p_start;		// starting point of the frame as index of the image
} framed_t;

static inline int min(const int a, const int b) {return((a < b) ? a : b);}
static inline int max(const int a, const int b) {return((a > b) ? a : b);}

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
			--i;
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
	//DEBUG_VAR(chunk_w);
	//DEBUG_VAR(chunk_h);
	
	int thed_w = t_w - w + chunk_w * t_w - 1;
	int thed_h = t_h - h + chunk_h * t_h - 1;
	//DEBUG_VAR(thed_w);
	//DEBUG_VAR(thed_h);

	int sw = chunk_w + (j > thed_w);
	int sh = chunk_h + (i > thed_h);
	//DEBUG_VAR(sw);
	//DEBUG_VAR(sh);
	
	int tl_j = chunk_w * j + max(0, j - thed_w - 1);
	int tl_i = chunk_h * i + max(0, i - thed_h - 1);
	int br_j = tl_j + sw - 1;
	int br_i = tl_i + sh - 1;
	//DEBUG_VAR(tl_j);
	//DEBUG_VAR(tl_i);
	//DEBUG_VAR(br_j);
	//DEBUG_VAR(br_i);

	int offset_l = min(halo, tl_j);
	int offset_r = min(halo, w - br_j - 1);
	int offset_t = min(halo, tl_i);
	int offset_b = min(halo, h - br_i - 1);
	//DEBUG_VAR(offset_l);
	//DEBUG_VAR(offset_r);
	//DEBUG_VAR(offset_t);
	//DEBUG_VAR(offset_b);
	
	sw += offset_l + offset_r;
	sh += offset_t + offset_b;
	
	int p_start = (tl_j - offset_l) + w * (tl_i - offset_t);

	frame->sw = sw;
	frame->sh = sh;
	frame->p_start = p_start;
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

	image_t img = {.w = w, .h = h};

	// Here we want to implement the touch-by-one policy.
	// However in order to simulate a real situation when
	// data arrive from outside, even in the touch-by-all
	// policy the matrix is initialized by the master thread
	// then one has to implement a local copy policy. Still
	// we have memory contention but this time only in the
	// initializing phase. The real computation is free of
	// this issue. However it requires the double of memory.
	// The question in: all in all what is the best way?
	// This is just an experiment to simulate the image blurring.
	
	uint16_t* data = (uint16_t*)calloc(w * h, sizeof(uint16_t));
	if (!data) {
		fprintf(stderr, "Bad allocation: not enough memory.\n");
		return 1;
	}

	{
	        long pseudo_seed = (long)(strlen(getenv("PATH")) * getpid())^(long)(&argc);
		struct drand48_data rnd;
		double pseudo_random;
	        srand48_r(pseudo_seed, &rnd);

		for (int i = 0; i < w * h; ++i)
			data[i] = (drand48_r(&rnd, &pseudo_random),(uint16_t)(pseudo_random * UINT16_MAX + 0.5));

		printf("Data initialized by master.\n");
	}

	int threads;
	thread_image_t t_img;

	
	#pragma omp parallel
	{
		#pragma omp master
		{
			threads = omp_get_num_threads();
			printf("threads number: %d\n", threads);
			printf("allocated memory: %f Mb.\n", ((float)(w*h*sizeof(int16_t))) / ((float)(2 << 19)));
			grid_dimension(&img, threads, &t_img);
			printf("w_threads: %d, h_threads: %d.\n\n", t_img.t_w, t_img.t_h);	
		}

		#pragma omp barrier

		
		int myid = omp_get_thread_num();
		int i = myid / t_img.t_w;
		int j = myid % t_img.t_w;
		framed_t frame;
		//get_frame(&t_img, k, i, j, &frame);

		#pragma omp critical(show)
		{
			get_frame(&t_img, k, i, j, &frame);
			printf("[THREAD %d] coordinates: (%d, %d)\n"
			       "            dimension: (%d, %d)\n"
			       "            p_start: %d\n",
			       myid, i, j, frame.sw, frame.sh, frame.p_start);
		}
	     	
	}
	
	free(data);
	
	
	return 0;
}
