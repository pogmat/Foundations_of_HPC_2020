#include <malloc.h>
#include "kernel.h"

kernel_t*  init_kernel_from_vector(const real* const init_vector,
			    const unsigned int s)
{
	// Warning it assumes that vectro lenght should be
	// exactly (2s+1)^2.

	kernel_t* new_kernel =(kernel_t*)malloc(sizeof(kernel_t));
	if (!new_kernel) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}

	new_kernel->s = s;
	
	unsigned int kernel_size = (2*s + 1) * (2*s + 1); 
	new_kernel->kernel = (real*)calloc(kernel_size, sizeof(real));
	if (!new_kernel->kernel) {
		fprintf(stderr, "Bad kernel allocaton.\n");
		return NULL;
	}

	if (init_vector) {
		for (unsigned int i = 0; i < kernel_size; ++i) {
			new_kernel->kernel[i] = init_vector[i];
		}
	}
	
	return new_kernel;
}

kernel_t* init_kernel_from_file(const char* const filename)
{
	// The first number of the file must be s.
	// It mus be followed by (2s+1)^2 numbers.

	kernel_t* new_kernel =(kernel_t*)malloc(sizeof(kernel_t));
	if (!new_kernel) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}
	
	FILE* fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "No such file: %s\n", filename);
		return NULL;
	}

	unsigned int s;
	
	if (fscanf(fp, "%u", &s) != 1) {
		fprintf(stderr, "Bad kernel file format: wrong first line.\n");
		fclose(fp);
		return NULL;
	}

	new_kernel->s = s;
	
	unsigned int kernel_size = (2*s + 1) * (2*s + 1);
	new_kernel->kernel = (real*)calloc(kernel_size, sizeof(real));
	if (!new_kernel) {
		fprintf(stderr, "Bat allocation.\n");
		fclose(fp);
		return NULL;
	}
	
	unsigned int i = 0;
	int read;
	real buffer;
	// TODO: It is not portable because if real is float %lf -> %f
	while ((read = fscanf(fp, "%lf", &buffer)) != EOF) {
		if (!read) {
			fprintf(stderr, "Bad kernel file at line %d", i + 1);
			fclose(fp);
			free(new_kernel->kernel);
			return NULL;
		}
		new_kernel->kernel[i] = buffer;
		++i;
	}

	if (i != kernel_size) {
		fprintf(stderr, "Bad kernel file: wrong dimensions.\n");
		fclose(fp);
		return NULL;
	}

	fclose(fp);
	
	return new_kernel;
}

kernel_t* copy_kernel(const kernel_t* const k)
{
	kernel_t* new_kernel =(kernel_t*)malloc(sizeof(kernel_t));
	if (!new_kernel) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}

	new_kernel->kernel = (real*)malloc(k->s * sizeof(real));
	if (!new_kernel->kernel) {
		fprintf(stderr, "Bad kernel allocaton.\n");
		return NULL;
	}

	for (unsigned int i = 0; i < k->s; ++i)
		new_kernel->kernel[i] = k->kernel[i];

	new_kernel->s = k->s;
	
	return new_kernel;
}

real get_luminosity(const kernel_t* const k)
{
	real luminosity = 0;
	unsigned int kernel_size = (2 * k->s + 1) * (2 * k->s +1);
	for (unsigned int i = 0; i < kernel_size; ++i)
		luminosity += k->kernel[i];

	return luminosity;
}

real get_partial_luminosity(const kernel_t* const k,
			    const int top,
			    const int bottom,
			    const int left,
			    const int right)
{
	real luminosity = 0;
	int kernel_width = 2 * k->s + 1;
	for (int i = top + (int)k->s; i <= bottom + (int)k->s; ++i)
		for (int j = left + (int)k->s; j <= right + (int)k->s; ++j)
			luminosity += k->kernel[j + kernel_width * i];
       
	return luminosity;
}

void normalize_luminosity(kernel_t* const k, const real norm)
{
	real luminosity = get_luminosity(k);
	real ratio = norm / luminosity;
	unsigned int kernel_size = (2 * k->s + 1) * (2 * k->s +1);
	for (unsigned int i = 0; i < kernel_size; ++i)
		k->kernel[i] *= ratio;
}

void free_kernel(kernel_t** const k)
{
	if (!(*k))
		return;
	
	free((*k)->kernel);
	free(*k);
	*k = NULL;
}
