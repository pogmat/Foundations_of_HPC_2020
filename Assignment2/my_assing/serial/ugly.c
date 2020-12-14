#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "../pgm/pgm_bin.h"
#include "../kernel/kernel.h"

static inline int max(const int a, const int b) {return ((a < b) ? b : a);}
static inline int min(const int a, const int b) {return ((a < b) ? a : b);}

int blur_pgm(const pgm_file* const original,
	     const kernel_t* const kernel,
	     pgm_file* const new)
{
	strcpy(new->magic, original->magic);
	register int w = new->width = original->width;
	register int h = new->height = original->height;
	new->maximum_value = original->maximum_value;

        size_t bytes_data = ((original->maximum_value > UINT8_MAX) ? 2 : 1)
		* original->width * original->height;

	new->data = (byte*)malloc(bytes_data * sizeof(byte));
	if (!new->data) {
		fprintf(stderr, "Bad allocation.\n");
		return -1;
	}

	register int s = kernel->s;
	register real new_value;
	
	if (original->maximum_value > UINT8_MAX) {
		dbyte* n_p = (dbyte*)new->data;
		dbyte* o_p = (dbyte*)original->data;
		real* k_p = kernel->kernel;

		for (int i = 0; i < h; ++i)
			for (int j = 0; j < w; ++j) {
				new_value = 0.0;
				for (int a = max(0, i - s); a < min(h, i + s + 1); ++a)
					for (int b = max(0, j - s); b < min(w, j + s +1); ++b)
						new_value +=
							o_p[b + w * a] * k_p[b - j + s + (2 * s +1) * (a - i + s)];
				n_p[j + w * i] = (dbyte)min(UINT16_MAX, (uint64_t)(new_value + 0.5));
			}		

	} else {
		byte* n_p = new->data;
		byte* o_p = original->data;
		real* k_p = kernel->kernel;

		for (int i = 0; i < h; ++i)
			for (int j = 0; j < w; ++j) {
				new_value = 0.0;
				for (int a = max(0, i - s); a < min(h, i + s + 1); ++a)
					for (int b = max(0, j - s); b < min(w, j + s +1); ++b)
						new_value +=
							o_p[b + w * a] * k_p[b - j + s + (2 * s +1) * (a - i + s)];
				n_p[j + w * i] = (byte)min(UINT8_MAX, (uint64_t)(new_value + 0.5));
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
