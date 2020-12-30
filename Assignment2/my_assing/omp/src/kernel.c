#include <malloc.h>
#include <math.h>
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

kernel_t* mean_kernel(const unsigned int s)
{
	kernel_t* new_kernel =(kernel_t*)malloc(sizeof(kernel_t));
	if (!new_kernel) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}

	new_kernel->s = s;
	
	const unsigned int kernel_size = (2*s + 1) * (2*s + 1);
	new_kernel->kernel = (real*)calloc(kernel_size, sizeof(real));
	if (!new_kernel->kernel) {
		fprintf(stderr, "Bad kernel allocaton.\n");
		return NULL;
	}
	
	const double val = 1 / ((double)kernel_size);
	for (unsigned int i = 0; i < kernel_size; ++i)
		new_kernel->kernel[i] = val;
	
	return new_kernel;
}

kernel_t* weight_kernel(const unsigned int s, const double f)
{
	kernel_t* new_kernel =(kernel_t*)malloc(sizeof(kernel_t));
	if (!new_kernel) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}

	new_kernel->s = s;
	
	const unsigned int kernel_size = (2*s + 1) * (2*s + 1);
	new_kernel->kernel = (real*)calloc(kernel_size, sizeof(real));
	if (!new_kernel->kernel) {
		fprintf(stderr, "Bad kernel allocaton.\n");
		return NULL;
	}	
	
	const double w = ((double)(1 - f)) / ((double)(kernel_size - 1));
	for (unsigned int i = 0; i < kernel_size / 2; ++i)
		new_kernel->kernel[i] = w;
	new_kernel->kernel[kernel_size / 2] = f;
	for (unsigned int i = kernel_size / 2 + 1; i < kernel_size; ++i)
		new_kernel->kernel[i] = w;
		
	return new_kernel;
}

kernel_t* gaussian_kernel(const unsigned int s)
{
	kernel_t* new_kernel =(kernel_t*)malloc(sizeof(kernel_t));
	if (!new_kernel) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}

	new_kernel->s = s;
	
	const unsigned int kernel_size = (2*s + 1) * (2*s + 1);
	new_kernel->kernel = (real*)calloc(kernel_size, sizeof(real));
	if (!new_kernel->kernel) {
		fprintf(stderr, "Bad kernel allocaton.\n");
		return NULL;
	}
	
	const double invs2 = 1 / ((double)(s * s));
	for (int i = -s; i <= (int)s; ++i)
		for (int j = -s; j <= (int)s; ++j)
			new_kernel->kernel[j + s + (2*s+1) * (i + s)] = exp(-((double)(i*i + j*j)) * invs2);

	normalize_luminosity(new_kernel, 1.0);
	
	return new_kernel;
}

kernel_t* copy_kernel(const kernel_t* const k)
{
	kernel_t* new_kernel =(kernel_t*)malloc(sizeof(kernel_t));
	if (!new_kernel) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}

	new_kernel->kernel = (real*)malloc((2 * k->s + 1) * (2 * k->s + 1) * sizeof(real));
	if (!new_kernel->kernel) {
		fprintf(stderr, "Bad kernel allocaton.\n");
		return NULL;
	}

	for (unsigned int i = 0; i < (1 + 2 * k->s) * (1 + 2 * k->s); ++i)
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
