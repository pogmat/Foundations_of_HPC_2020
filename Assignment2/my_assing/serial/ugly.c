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

int remove_vignetting_inplace(const pgm_file* const original,
			      const kernel_t* const kernel)
// WARNING: If we want to do the implace trasformation without allocate further memory
//	    we cannot declare const byte* data. Then without this const it is of course
//          possible to modify data even if the pgm_file is pass as const.


// THINK:   Perhaps it would be better to normalize every pixel directly inside blur.
{
	// NOTE: we split the boder in 4 stripes + 4 corners:
	//       for the stripes the renormalization factor has to be computed only once,
	//       for corners has to be computed for every pixel.

	int w = original->width;
	int h = original->height;
	int halo = kernel->s;
	real total_luminosity = get_luminosity(kernel);
	real luminosity;

        if (original->maximum_value > UINT8_MAX) {
		dbyte* o_p = (dbyte*)original->data;
		
		// Top
		for (int i = halo - 1; i >= 0; --i) {
			luminosity = get_partial_luminosity(kernel, max(-halo, -i), min(halo, h - i - 1), -halo, halo);
			if (!luminosity)
				return -1;
			for (int j = halo; j < w - halo; ++j) {
				o_p[i * w + j] = (uint16_t)min(UINT16_MAX,
							       (o_p[i * w + j] *  total_luminosity / luminosity + 0.5));
			}
		}
		
		// Bottom
		// the max here is to prevent overlap with the Top vignetting removal
		for (int i = max(h - halo, halo); i < h; ++i) {
			luminosity = get_partial_luminosity(kernel, max(-halo, -i), min(halo, h - i - 1), -halo, halo);
			if (!luminosity)
				return -1;
			for (int j = halo; j < w - halo; ++j)
				o_p[i * w + j] = (uint16_t)min(UINT16_MAX,
							       (o_p[i * w + j] *  total_luminosity / luminosity + 0.5));
		}

		// Left
		for (int j = halo - 1; j >= 0; --j) {
			luminosity = get_partial_luminosity(kernel, -halo, halo, max(-halo, -j), min(halo, w - j - 1));
			if (!luminosity)
				return -1;
			for (int i = halo; i < h - halo; ++i)
				o_p[i * w + j] = (uint16_t)min(UINT16_MAX,
							       (o_p[i * w + j] *  total_luminosity / luminosity + 0.5));
		}

		// Right
		// the max here is to prevent overlap with the Bottom vignetting removal
		for (int j = max(w - halo, halo); j < w; ++j) {
			luminosity = get_partial_luminosity(kernel, -halo, halo, max(-halo, -j), min(halo, w - j - 1));
			if (!luminosity)
				return -1;
			for (int i = halo; i < h - halo; ++i)
				o_p[i * w + j] = (uint16_t)min(UINT16_MAX,
						     (o_p[i * w + j] *  total_luminosity / luminosity + 0.5));
	       	}
	     
		// Corners
		int top, bottom;
		int i_second_min = max(h - halo, halo);	// to avoid overlap with other corners
		int j_second_min = max(w - halo, halo); // to avoid overlap with other corners
		for (int i = 0; i < halo; ++i) {
			top = max(-halo, -i);
			bottom = min(halo, h - i - 1);
			for (int j = 0; j < halo; ++j) {
				luminosity =  get_partial_luminosity(kernel,
								     top,
								     bottom,
								     max(-halo, -j),
								     min(halo, w - j - 1));
				if (!luminosity)
					return -1;
				o_p[i * w + j] = (uint16_t)min(UINT16_MAX,
							       (o_p[i * w + j] * total_luminosity / luminosity + 0.5));
			}
			for (int j = j_second_min; j < w; ++j) {
				luminosity =  get_partial_luminosity(kernel,
								     top,
								     bottom,
								     max(-halo, -j),
								     min(halo, w - j - 1));
				if (!luminosity)
					return -1;
				o_p[i * w + j] = (uint16_t)min(UINT16_MAX,
							       (o_p[i * w + j] * total_luminosity / luminosity + 0.5));
			}
		     		}
		for (int i = i_second_min; i < h; ++i) {
			top = max(-halo, -i);
			bottom = min(halo, h - i - 1);
			for (int j = 0; j < halo; ++j) {
				luminosity =  get_partial_luminosity(kernel,
								     top,
								     bottom,
								     max(-halo, -j),
								     min(halo, w - j - 1));
				if (!luminosity)
					return -1;
				o_p[i * w + j] = (uint16_t)min(UINT16_MAX,
							       (o_p[i * w + j] * total_luminosity / luminosity + 0.5));
			}
			for (int j = j_second_min; j < w; ++j) {
				luminosity =  get_partial_luminosity(kernel,
								     top,
								     bottom,
								     max(-halo, -j),
								     min(halo, w - j - 1));
				if (!luminosity)
					return -1;
				o_p[i * w + j] = (uint16_t)min(UINT16_MAX,
							       (o_p[i * w + j] * total_luminosity / luminosity + 0.5));
			}
		}
	} else {
		byte* o_p = original->data;

		// Top
		for (int i = halo - 1; i >= 0; --i) {
			luminosity = get_partial_luminosity(kernel, max(-halo, -i), min(halo, h - i - 1), -halo, halo);
			if (!luminosity)
				return -1;
			for (int j = halo; j < w - halo; ++j) {
				o_p[i * w + j] = (uint8_t)min(UINT8_MAX,
							      (o_p[i * w + j] *  total_luminosity / luminosity + 0.5));
			}
		}

		// Bottom
		// the max here is to prevent overlap with the Top vignetting removal
		for (int i = max(h - halo, halo); i < h; ++i) {
			luminosity = get_partial_luminosity(kernel, max(-halo, -i), min(halo, h - i - 1), -halo, halo);
			if (!luminosity)
				return -1;
			for (int j = halo; j < w - halo; ++j)
				o_p[i * w + j] = (uint8_t)min(UINT8_MAX,
							      (o_p[i * w + j] *  total_luminosity / luminosity + 0.5));
		}

		// Left
		for (int j = halo - 1; j >= 0; --j) {
			luminosity = get_partial_luminosity(kernel, -halo, halo, max(-halo, -j), min(halo, w - j - 1));
			if (!luminosity)
				return -1;
			for (int i = halo; i < h - halo; ++i)
				o_p[i * w + j] = (uint8_t)min(UINT8_MAX,
							      (o_p[i * w + j] *  total_luminosity / luminosity + 0.5));
		}

		// Right
		// the max here is to prevent overlap with the Bottom vignetting removal
		for (int j = max(w - halo, halo); j < w; ++j) {
			luminosity = get_partial_luminosity(kernel, -halo, halo, max(-halo, -j), min(halo, w - j - 1));
			if (!luminosity)
				return -1;
			for (int i = halo; i < h - halo; ++i)
				o_p[i * w + j] = (uint8_t)min(UINT8_MAX,
							      (o_p[i * w + j] *  total_luminosity / luminosity + 0.5));
	       	}

		// Corners
		int top, bottom;
		int i_second_min = max(h - halo, halo);	// to avoid overlap with other corners
		int j_second_min = max(w - halo, halo); // to avoid overlap with other corners
		for (int i = 0; i < halo; ++i) {
			top = max(-halo, -i);
			bottom = min(halo, h - i - 1);
			for (int j = 0; j < halo; ++j) {
				luminosity =  get_partial_luminosity(kernel,
								     top,
								     bottom,
								     max(-halo, -j),
								     min(halo, w - j - 1));
				if (!luminosity)
					return -1;
				o_p[i * w + j] = min(UINT8_MAX,
						     (uint8_t)(o_p[i * w + j] * total_luminosity / luminosity + 0.5));
			}
			for (int j = j_second_min; j < w; ++j) {
				luminosity =  get_partial_luminosity(kernel,
								     top,
								     bottom,
								     max(-halo, -j),
								     min(halo, w - j - 1));
				if (!luminosity)
					return -1;
				o_p[i * w + j] = (uint8_t)min(UINT8_MAX,
							      (o_p[i * w + j] * total_luminosity / luminosity + 0.5));
			}
		}
		for (int i = i_second_min; i < h; ++i) {
			top = max(-halo, -i);
			bottom = min(halo, h - i - 1);
			for (int j = 0; j < halo; ++j) {
				luminosity =  get_partial_luminosity(kernel,
								     top,
								     bottom,
								     max(-halo, -j),
								     min(halo, w - j - 1));
				if (!luminosity)
					return -1;
				o_p[i * w + j] = (uint8_t)min(UINT8_MAX,
							      (o_p[i * w + j] * total_luminosity / luminosity + 0.5));
			}
			for (int j = j_second_min; j < w; ++j) {
				luminosity =  get_partial_luminosity(kernel,
								     top,
								     bottom,
								     max(-halo, -j),
								     min(halo, w - j - 1));
				if (!luminosity)
					return -1;
				o_p[i * w + j] = (uint8_t)min(UINT8_MAX,
							      (o_p[i * w + j] * total_luminosity / luminosity + 0.5));
			}
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

	#ifndef KEEP_VIGNETTING
	if (remove_vignetting_inplace(&new_pgm, &my_kernel)) {
		fprintf(stderr, "Problem with vignetting.\n");
		free_kernel(&my_kernel);
		free_pgm(&original_pgm);
		free_pgm(&new_pgm);		
	}
	#endif
	
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
