#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "pgm_bin.h"
#include "kernel.h"
#include "grid.h"

#define _GNU_SOURCE

#ifndef __GNUC__
#error This program require GNU C extensions
#endif
#define CLEANUP(x) __attribute__((__cleanup__(x)))

#define TIMESTAMP_START t_start
#define TIMESTAMP_STOP t_stop
#define T_START clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(TIMESTAMP_START))
#define T_STOP clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(TIMESTAMP_STOP))
#define DELTA_mSEC 1e3 * (double)(TIMESTAMP_STOP.tv_sec - TIMESTAMP_START.tv_sec) + 1e-6 * (double)(TIMESTAMP_STOP.tv_nsec - TIMESTAMP_START.tv_nsec) 

#ifndef _OPENMP
#warning You are compiling without openmp support!
#endif

void cleanup_pgm_file(pgm_file_t** f) {close_pgm(f);}
void cleanup_kernel_t(kernel_t** k) {free_kernel(k);}
