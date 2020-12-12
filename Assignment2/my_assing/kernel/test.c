#include <stdio.h>
#include "../pgm/pgm_bin.h"
#include "../kernel/kernel.h"

int main(int argc, char** argv)
{
	if (argc != 2) {
		fprintf(stderr, "Give me a kernel file.\n");
		return 1;
	}
	
	kernel_t my_kernel;

	if (init_kernel_from_file(argv[1], &my_kernel)) {
		fprintf(stderr, "Problem with kernel file reading.\n");
		return 1;
	}

	real luminosity = get_luminosity(&my_kernel);
	printf("luminosity = %lf\n", luminosity);
	normalize_luminosity(&my_kernel, 1.);

	printf("s: %d\n", my_kernel.s);
	for (int i = 0; i < 9; ++i)
		printf("kernel[%d] = %f\n", i, my_kernel.kernel[i]);

	free_kernel(&my_kernel);
	
	return 0;
}
