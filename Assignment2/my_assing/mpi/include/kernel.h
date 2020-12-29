#ifndef _KERNEL_H
#define _KERNEL_H

typedef double real;

typedef struct {
	unsigned int s;
	real* kernel;
} kernel_t;


kernel_t* init_kernel_from_vector(const real* const, const unsigned int);
kernel_t* init_kernel_from_file(const char* const);
kernel_t* mean_kernel(const unsigned int);
kernel_t* weight_kernel(const unsigned int, const double);
kernel_t* gaussian_kernel(const unsigned int);
kernel_t* copy_kernel(const kernel_t* const);
real get_luminosity(const kernel_t* const);
real get_partial_luminosity(const kernel_t* const, const int, const int, const int, const int);
void normalize_luminosity(kernel_t* const, const real);
void free_kernel(kernel_t** const);

#endif
