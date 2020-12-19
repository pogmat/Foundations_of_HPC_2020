#include <stdio.h>
#include <stdlib.h>
#include "../pgm/pgm_bin.h"
#include "../kernel/kernel.h"

int main(int argc, char** argv)
{
	if (argc != 6) {
		fprintf(stderr, "Give me a kernel file and four integer.\n");
		return 1;
	}

	int top = atoi(argv[2]);
	int bottom = atoi(argv[3]);
	int left = atoi(argv[4]);
	int right = atoi(argv[5]);
	
	kernel_t my_kernel;

	if (init_kernel_from_file(argv[1], &my_kernel)) {
		fprintf(stderr, "Problem with kernel file reading.\n");
		return 1;
	}

	if ((abs(top) > my_kernel.s) ||
	    (abs(bottom) > my_kernel.s) ||
	    (abs(left) > my_kernel.s) ||
	    (abs(right) > my_kernel.s)) {
		    fprintf(stderr, "Wrong cut dimension.\n");
		    free_kernel(&my_kernel);
		    return 1;
	    }
	
	real luminosity = get_luminosity(&my_kernel);
	printf("total luminosity = %lf\n", luminosity);
	normalize_luminosity(&my_kernel, 1.);
	real partial_luminosity = get_partial_luminosity(&my_kernel, top, bottom, left, right);
		printf("partial luminosity = %lf\n", partial_luminosity);

	printf("s: %d\n", my_kernel.s);
	for (int i = 0; i < 9; ++i)
		printf("kernel[%d] = %f\n", i, my_kernel.kernel[i]);

	free_kernel(&my_kernel);
	
	return 0;
}
