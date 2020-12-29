#include <stdio.h>
#include <stdlib.h>
#include "kernel.h"

#define CLEANUP(x) __attribute__((__cleanup__(x)))

void cleanup_kernel_t(kernel_t** kp){free_kernel(kp);}

void print_kernel(const kernel_t* k)
{
	const int s = k->s;

	for (int i = -s; i <= (int)s; ++i) {
		for (int j = -s; j <= (int)s; ++j)
			printf("%7.6lf ", k->kernel[j + s + (2*s+1) * (i + s)]);
		printf("\n");
	}
			 
}

int main(/*int argc, char** argv*/)
{
	/*
	if (argc != 6) {
		fprintf(stderr, "Give me a kernel file and four integer.\n");
		return 1;
	}

	int top = atoi(argv[2]);
	int bottom = atoi(argv[3]);
	int left = atoi(argv[4]);
	int right = atoi(argv[5]);
	
	kernel_t* CLEANUP(cleanup_kernel_t) kernel1 = NULL;
	kernel_t* CLEANUP(cleanup_kernel_t) kernel2 = NULL;
	
	kernel1 = init_kernel_from_file(argv[1]);
	if (!kernel1) {
		fprintf(stderr, "Problem with kernel file reading.\n");
		return 1;
	}

	if ((abs(top) > kernel1->s) ||
	    (abs(bottom) > kernel1->s) ||
	    (abs(left) > kernel1->s) ||
	    (abs(right) > kernel1->s)) {
		    fprintf(stderr, "Wrong cut dimension.\n");
		    return 1;
	    }
	
	real luminosity = get_luminosity(kernel1);
	printf("total luminosity = %lf\n", luminosity);
	normalize_luminosity(kernel1, 1.);
	real partial_luminosity = get_partial_luminosity(kernel1, top, bottom, left, right);
	printf("partial luminosity = %lf\n", partial_luminosity);
		
	printf("s: %d\n", kernel1->s);
	for (unsigned int i = 0; i < (1+2*kernel1->s)*(1+2*kernel1->s); ++i)
		printf("kernel[%d] = %f\n", i, kernel1->kernel[i]);

	kernel2 = copy_kernel(kernel1);
	if (!kernel2) {
		fprintf(stderr, "Problem with kernel file reading.\n");
		return 1;
	}	
	printf("s: %d\n", kernel2->s);
	for (unsigned int i = 0; i < kernel1->s; ++i)
		printf("kernel[%d] = %f\n", i, kernel2->kernel[i]);
	*/

	kernel_t* CLEANUP(cleanup_kernel_t) m = mean_kernel(5);;
	kernel_t* CLEANUP(cleanup_kernel_t) w = weight_kernel(5, 0.52);
	kernel_t* CLEANUP(cleanup_kernel_t) g = gaussian_kernel(5);	

	print_kernel(m);
	printf("\n");
	print_kernel(w);
	printf("\n");
	print_kernel(g);
	
	return 0;
}
