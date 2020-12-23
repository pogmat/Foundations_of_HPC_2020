#ifndef _KERNEL_H
#define _KERNEL_H

typedef double real;

typedef struct {
	unsigned int s;
	real* kernel;
} kernel_t;

int init_kernel_from_vector(const real* const, const unsigned int, kernel_t* const);
int init_kernel_from_file(const char* const, kernel_t* const);
real get_luminosity(const kernel_t* const);
real get_partial_luminosity(const kernel_t* const, const int, const int, const int, const int);
void normalize_luminosity(kernel_t* const, const real);
void free_kernel(kernel_t* const);

#endif
