#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "../pgm/pgm_bin.h"
#include "../kernel/kernel.h"

int blur_pgm(const pgm_file* const original,
	     const kernel_t* const kernel,
	     pgm_file* const new)
{
	strcpy(new->magic, original->magic);
	register int w = new->width = original->width;
	register int h = new->height = original->height;
	new->maximum_value = original->maximum_value;

        size_t bytes_data = ((original->maximum_value > 255) ? 2 : 1)
		* original->width * original->height;

	new->data = (byte*)malloc(bytes_data * sizeof(byte));
	if (!new->data) {
		fprintf(stderr, "Bad allocation.\n");
		return -1;
	}

	register int s = kernel->s;
	register int new_value;
	
	if (original->maximum_value > 255) {
		unsigned short* new_p = (unsigned short*)new->data;
		unsigned short* original_p = (unsigned short*)original->data;

		for (int i = 0; i < w; ++i)
			for (int j = 0; j < h; ++j) {
				new_value = 0;
				for (int a = -s; a <= s; ++a)
					for (int b = -s; b <= s; ++b) {
						if (((i + a) >= 0) &&
						    ((i + a) <  w) &&
						    ((j + b) >= 0) &&
						    ((j + b) <  h)) {
							new_value +=(unsigned int)
								(kernel->kernel[(a + s) + (b + s) * (2 * s + 1)]
								* original_p[(i + a) + (j + b) * w] + 0.5);
						}
					}
				new_p[i + j * w] = (unsigned short)(new_value > 255 ? 255 : new_value);
			}	       
	} else {
		for (int i = 0; i < w; ++i)
			for (int j = 0; j < h; ++j) {
				new_value = 0;
				for (int a = -s; a <= s; ++a)
					for (int b = -s; b <= s; ++b) {
						if (((i + a) >= 0) &&
						    ((i + a) <  w) &&
						    ((j + b) >= 0) &&
						    ((j + b) <  h)) {
							new_value +=(unsigned int)
								(kernel->kernel[(a + s) + (b + s) * (2 * s + 1)]
								* original->data[(i + a) + (j + b) * w] + 0.5);
						}
					}
				new->data[i + j * w] = (byte)(new_value > 255 ? 255 : new_value);
				}
	}
	
	return 0;
}

int main(int argc, char** argv)
{
	if (argc != 4) {
		fprintf(stderr,
			"USEGE:\n%s "
		        "<input_image> "
			"<input_kernel> "
			"<output_image>\n",
			argv[0]);

		return 1;
	}
	
	kernel_t my_kernel;
	if (init_kernel_from_file(argv[2], &my_kernel)) {
		fprintf(stderr, "Problem with kernel file reading.\n");
		return 1;
	}
	normalize_luminosity(&my_kernel, 1.0);

	pgm_file original_pgm, new_pgm;
	if (read_pgm(argv[1], &original_pgm)) {
		fprintf(stderr, "Problem with pgm file reading.\n");
		free_kernel(&my_kernel);
		return 1;
	}

	if (blur_pgm(&original_pgm, &my_kernel, &new_pgm)) {
		fprintf(stderr, "Problem with blurring procedure.\n");
		free_kernel(&my_kernel);
		free_pgm(&original_pgm);
		return 1;
	}

	if (write_pgm(argv[3], &new_pgm)) {
		fprintf(stderr, "Problem with pgm file writing.\n");
		free_kernel(&my_kernel);
		free_pgm(&original_pgm);
		free_pgm(&new_pgm);
	}

	free_kernel(&my_kernel);
	free_pgm(&original_pgm);
	free_pgm(&new_pgm);
	
	return 0;
}
