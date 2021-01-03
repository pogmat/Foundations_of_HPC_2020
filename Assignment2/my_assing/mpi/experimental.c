#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>


#ifndef __GNUC__
#error This uses GNU C extensions
#endif

#define MASTER 0

#define CLEANUP(x) __attribute__((__cleanup__(x)))

typedef struct {
	char magic[3];
	int w;
	int h;
	int dept;
	long offset;
} metadata_t;

typedef struct {
	int w;
	int h;
	int j;
	int i;
} tharea_t;

typedef struct {
	tharea_t read;
	tharea_t writ;
} thimg_t;

typedef struct {
	unsigned int s;
	double* kernel;
} kernel_t;

static inline int min(const int a, const int b) {return((a < b) ? a : b);}
static inline int max(const int a, const int b) {return((a > b) ? a : b);}

double get_luminosity(const kernel_t* const k)
{
	double luminosity = 0;
	unsigned int kernel_size = (2 * k->s + 1) * (2 * k->s +1);
	for (unsigned int i = 0; i < kernel_size; ++i)
		luminosity += k->kernel[i];

	return luminosity;
}

double get_partial_luminosity(const kernel_t* const k,
			    const int top,
			    const int bottom,
			    const int left,
			    const int right)
{
	double luminosity = 0;
	int kernel_width = 2 * k->s + 1;
	for (int i = top + (int)k->s; i <= bottom + (int)k->s; ++i)
		for (int j = left + (int)k->s; j <= right + (int)k->s; ++j)
			luminosity += k->kernel[j + kernel_width * i];
       
	return luminosity;
}

void normalize_luminosity(kernel_t* const k, const double norm)
{
	double luminosity = get_luminosity(k);
	double ratio = norm / luminosity;
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

kernel_t* mean_kernel(const unsigned int s)
{
	kernel_t* new_kernel =(kernel_t*)malloc(sizeof(kernel_t));
	if (!new_kernel) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}

	new_kernel->s = s;
	
	const unsigned int kernel_size = (2*s + 1) * (2*s + 1);
	new_kernel->kernel = (double*)calloc(kernel_size, sizeof(double));
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
	new_kernel->kernel = (double*)calloc(kernel_size, sizeof(double));
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
	new_kernel->kernel = (double*)calloc(kernel_size, sizeof(double));
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

void cleanup_kernel_t(kernel_t** k) {free_kernel(k);}
void cleanup_file(FILE** fpp) {fclose(*fpp);}
void cleanup_uint16_t(uint16_t** pp)
{
	if (*pp)
		free(*pp);

	*pp = NULL;
}

static void skip_comment(FILE * const fp)
{
	int ch;
	int inside_comment = 0;
	
	while ((ch = fgetc(fp)) != EOF) {
		switch (ch) {
		case ' ':
		case '\t':
			break;
		case '#':
			inside_comment = 1;
			break;
		case '\n':
			inside_comment = 0;
			break;
		default:
			if (!inside_comment) {
				fseek(fp, -1, SEEK_CUR);
				return;
			}
			break;
		}
	}
}


int pgm_get_metadata(const char* filename, metadata_t* md)
{
	FILE* CLEANUP(cleanup_file) fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "No such file: %s", filename);
		return -1;
	}

	skip_comment(fp);
	fgets(md->magic, sizeof(md->magic), fp);
	if (strncmp(md->magic, "P5", 3)) {
		fprintf(stderr, "Invalid format: only P5 format allowed.\n");
		return -1;
	}

	skip_comment(fp);
	if (!fscanf(fp, "%d", &md->w)) {
		fprintf(stderr, "Invalid format: no width.\n");
		return -1;
	}

	skip_comment(fp);
	if (!fscanf(fp, "%d", &md->h)) {
		fprintf(stderr, "Invalid format: no height.\n");
		return -1;
	}

	skip_comment(fp);
	if (!fscanf(fp, "%d", &md->dept)) {
		fprintf(stderr, "Invalid format: no dept.\n");
		return -1;
	}
	
	md->offset = ftell(fp) + 1;
	
	return 0;
}

int pgm_set_metadata(const char* filename, metadata_t* md)
{
	FILE* CLEANUP(cleanup_file) fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "No such file: %s", filename);
		return -1;
	}

	if (!fprintf(fp, "%2s\n%d %d\n%d\n",
		     md->magic,
		     md->w,
		     md->h,
		     md->dept)) {
		fprintf(stderr, "Impossible to write file header.\n");
		return -1;
	}
	
	md->offset = ftell(fp);
	
	return 0;
}


void grid_dimension(const int w, const int h,
		    const int threads,
		    int* const w_th, int* const h_th)
{
	int big_d, * big_p, * small_p;

	if (w > h) {
		big_d = w;
		big_p = w_th;
		small_p = h_th;
	} else {
		big_d = h;
		big_p = h_th;
		small_p = w_th;
	}

	int i = (int)((float)big_d / (sqrtf(w * h / threads)));

	while(1) {
		if (threads % i)
			++i;
		else {
			*big_p = i;
			*small_p = threads / i;
			break;
		}
	}
}

void get_frame(const int w, const int h,
	       const int w_th, const int h_th,
	       const int j, const int i,
	       const int halo,
	       thimg_t* const frame)
{
	int chunk_w = w / w_th;
	int chunk_h = h / h_th;

	int thed_w = w_th - w + chunk_w * w_th - 1;
	int thed_h = h_th - h + chunk_h * h_th - 1;

	frame->writ.w = chunk_w + (j > thed_w);
	frame->writ.h = chunk_h + (i > thed_h);

	int tl_j = chunk_w * j + max(0, j - thed_w - 1);
	int tl_i = chunk_h * i + max(0, i - thed_h - 1);
	int br_j = tl_j + frame->writ.w - 1;
	int br_i = tl_i + frame->writ.h - 1;

	frame->writ.j = tl_j;
	frame->writ.i = tl_i;

	int offset_l = min(halo, tl_j);
	int offset_r = min(halo, w - br_j - 1);
	int offset_t = min(halo, tl_i);
	int offset_b = min(halo, h - br_i - 1);

	frame->read.w = frame->writ.w + offset_l + offset_r;
	frame->read.h = frame->writ.h + offset_t + offset_b;
	frame->read.j = tl_j - offset_l;
	frame->read.i = tl_i - offset_t;
}

void endian_swap(uint16_t* data, size_t size)
{
	for (size_t i = 0; i < size; ++i)
		data[i] = ((data[i] & (uint16_t)0xff00) >> 8) + ((data[i] & (uint16_t)0x00ff) << 8);

}

uint16_t* blur_frame(const uint16_t* restrict original,
		     const thimg_t* restrict frame,
		     const kernel_t* restrict kernel)
{
	uint16_t* blurred = (uint16_t*)malloc(2 * frame->writ.w * frame->writ.h);
	if (!blurred) {
		fprintf(stderr, "Bad allocation.\n");
		return NULL;
	}
	double* k = kernel->kernel;
	register double t0, t1, t2, t3;
	register int i_offset, j_offset, o_offset, k_offset, b, b_start, b_stop;
	const int b_stride = 16;
	register int s = kernel->s;


	for (int i = 0; i < frame->writ.h; ++i)
		for (int j = 0; j < frame->writ.w; ++j) {
			t0 = t1 = t2 = t3 = 0;
			i_offset = frame->writ.i - frame->read.i;
			j_offset = frame->writ.j - frame->read.j;

			for (int a = max(0, i + i_offset - s); a < min(frame->read.h, i + i_offset + s + 1); ++a) {
				o_offset = frame->read.w * a;
				k_offset = - j - j_offset + s + (2 * s + 1) * (a - i - i_offset + s);
				b_start = max(0, j + j_offset - s);
				b_stop = min(frame->read.w, j + j_offset + s + 1);
				for (b = b_start; b < ((b_stop - b_start) / b_stride) * b_stride; b += b_stride) {
					t0 += (original[b +     o_offset] * k[b +     k_offset] +
					       original[b + 1 + o_offset] * k[b + 1 + k_offset] +
					       original[b + 2 + o_offset] * k[b + 2 + k_offset] +
					       original[b + 3 + o_offset] * k[b + 3 + k_offset]);
					t1 += (original[b + 4 + o_offset] * k[b + 4 + k_offset] +
					       original[b + 5 + o_offset] * k[b + 5 + k_offset] +
					       original[b + 6 + o_offset] * k[b + 6 + k_offset] +
					       original[b + 7 + o_offset] * k[b + 7 + k_offset]);
					t2 += (original[b + 8 + o_offset] * k[b + 8 + k_offset] +
					       original[b + 9 + o_offset] * k[b + 9 + k_offset] +
					       original[b +10 + o_offset] * k[b +10 + k_offset] +
					       original[b +11 + o_offset] * k[b +11 + k_offset]);					
					t3 += (original[b +12 + o_offset] * k[b +12 + k_offset] +
					       original[b +13 + o_offset] * k[b +13 + k_offset] +
					       original[b +14 + o_offset] * k[b +14 + k_offset] +
					       original[b +15 + o_offset] * k[b +15 + k_offset]);				       
				}

				for (; b < b_stop; ++b)
					t0 += original[b + o_offset] * k[b + k_offset];
			}
			blurred[j + frame->writ.w * i] = (uint16_t)((t0 + t1) + (t2 + t3));
		}

	return blurred;
}

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	int total, id;

	MPI_Comm_size(MPI_COMM_WORLD, &total);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	if (argc < 4) {
	usage:
		fprintf(stderr, "USAGE:\n%s [kernel-type] [kernel-size] {kernel-param} [input-file] {output-file}\n", argv[0]);
		return EXIT_FAILURE;
	}
	
	char default_out[] = "out.pgm";
	
	int kernel_type = atoi(argv[1]);
	unsigned int kernel_size = (unsigned int)atoi(argv[2]);
        char* output_filename = default_out;
	char* input_filename;
	double kernel_param;
	
	switch (kernel_type) {
	case 0:
	case 2:
		input_filename = argv[3];
		if (argc == 5)
			output_filename = argv[4];
		if (argc > 5)
			goto usage;
		break;
	case 1:
		if (argc < 5) {
			fprintf(stderr, "Weighted kernel needs one more parameter.\n");
			return EXIT_FAILURE;
		}
		kernel_param = atof(argv[3]);
		input_filename = argv[4];
		if (argc == 6)
			output_filename = argv[5];
		if (argc > 6)
			goto usage;
		break;
	default:
		fprintf(stderr, "Unkown kernel:\n1. Mean kernel\n2. Weighted kernel\n3. Gaussian kernel\n");
		return EXIT_FAILURE;
	}

	if (!(kernel_size % 2)) {
		fprintf(stderr, "Kernel size must be odd.\n");
		return EXIT_FAILURE;
	}

	int s = kernel_size / 2;

	kernel_t* CLEANUP(cleanup_kernel_t) kernel = NULL;

	switch (kernel_type) {
	case 0:
		kernel = mean_kernel((unsigned int)s);
		break;
	case 1:
		kernel = weight_kernel((unsigned int)s, kernel_param);
		break;
	case 2:
		kernel = gaussian_kernel((unsigned int)s);
	}

	if (!kernel) {
		fprintf(stderr, "Error in reading kernel file: %s.\n", argv[2]);
		return EXIT_FAILURE;
	}

	metadata_t metadata;
	MPI_Datatype MPI_metadata;
	{
		int len[5] = {3, 1, 1, 1, 1};
		MPI_Datatype types[5] = {MPI_CHAR, MPI_INT, MPI_INT, MPI_INT, MPI_LONG};
		MPI_Aint displ[5];
		MPI_Get_address(&metadata.magic, displ);
		MPI_Get_address(&metadata.w, displ + 1);
		MPI_Get_address(&metadata.h, displ + 2);
		MPI_Get_address(&metadata.dept, displ + 3);
		MPI_Get_address(&metadata.offset, displ + 4);
		MPI_Aint base = displ[0];
		for (int i = 0; i < 5; ++i)
			displ[i] = MPI_Aint_diff(displ[i], base);
		MPI_Type_create_struct(5, len, displ, types, &MPI_metadata);
		MPI_Type_commit(&MPI_metadata);
	}
	
	if (id == MASTER) {
		if(pgm_get_metadata(input_filename, &metadata)) {
			fprintf(stderr, "Problem in reading metadata in %s.\n", input_filename);
			EXIT_FAILURE;
		}
	}
	MPI_Bcast(&metadata, 1, MPI_metadata, MASTER, MPI_COMM_WORLD);
       	#ifdef DEBUG
	printf("[NODE %3d] w = %d, h = %d, dept = %d, offset = %ld\n",
	       id, metadata.w, metadata.h, metadata.dept, metadata.offset);
	#endif

	int th_axes[2];
	grid_dimension(metadata.w, metadata.h, total, th_axes, th_axes + 1);
	#ifdef DEBUG
	if (id == MASTER)
		printf("grid: %d x %d.\n", th_axes[0], th_axes[1]);
	#endif

	MPI_Comm CART_COMM;
	int periods[2] = {0, 0};
	MPI_Cart_create(MPI_COMM_WORLD, 2, th_axes, periods, 1, &CART_COMM);

	int coords[2];
	MPI_Cart_coords(CART_COMM, id, 2, coords);

	thimg_t frame;
	get_frame(metadata.w, metadata.h, th_axes[0], th_axes[1], coords[0], coords[1], kernel->s, &frame);

	#ifdef DEBUG
	printf("[NODE %3d] coord : (%d, %d)\n"
	       "           read d: (%d, %d)\n"
	       "           read p: (%d, %d)\n"
	       "           writ d: (%d, %d)\n"
	       "           writ p: (%d, %d)\n\n",
	       id, coords[0], coords[1],
	       frame.read.w, frame.read.h,
	       frame.read.j, frame.read.i,
	       frame.writ.w, frame.writ.h,
	       frame.writ.j, frame.writ.i);	
	#endif

	int array_dim[2] = {metadata.h, metadata.w};

	MPI_Datatype read_frame;
	{
		int subarray_dim[2] = {frame.read.h, frame.read.w};
		int start_point[2] = {frame.read.i, frame.read.j};
		MPI_Type_create_subarray(2,
					 array_dim,
					 subarray_dim,
					 start_point,
					 MPI_ORDER_C,
					 MPI_UINT16_T,
					 &read_frame);
		MPI_Type_commit(&read_frame);
	}

	//int pixel_size = 1 + metadata.dept > UINT8_MAX;
	uint16_t* in_data = (uint16_t*)malloc(2 * metadata.w * metadata.h);
	if (!in_data) {
		fprintf(stderr, "Bad allocation.\n");
		// TODO: do some cleanup
		MPI_Type_free(&read_frame);
		MPI_Type_free(&MPI_metadata);
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	MPI_File input_file;
	MPI_File_open(MPI_COMM_WORLD, input_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &input_file);
        MPI_File_set_view(input_file,
			  (MPI_Offset)metadata.offset,
			  MPI_UINT16_T,
			  read_frame,
			  "native",
			  MPI_INFO_NULL);
	MPI_File_read(input_file,
		      in_data,
		      frame.read.w * frame.read.h,
		      MPI_UINT16_T,
		      MPI_STATUS_IGNORE);
	endian_swap(in_data, frame.read.w * frame.read.h);
	
	
	if (id == MASTER) {
		if(pgm_set_metadata(output_filename, &metadata)) {
			fprintf(stderr, "Problem in writing metadata in %s.\n", input_filename);
			EXIT_FAILURE;
		}
	}
	MPI_Bcast(&metadata, 1, MPI_metadata, MASTER, MPI_COMM_WORLD);

	uint16_t* CLEANUP(cleanup_uint16_t) out_data = blur_frame(in_data, &frame, kernel); 
	
	MPI_Datatype writ_frame;
	{
		int subarray_dim[2] = {frame.writ.h, frame.writ.w};
		int start_point[2] = {frame.writ.i, frame.writ.j};
		MPI_Type_create_subarray(2,
					 array_dim,
					 subarray_dim,
					 start_point,
					 MPI_ORDER_C,
					  MPI_UINT16_T,
					 &writ_frame);
		MPI_Type_commit(&writ_frame);
	}
		
	
	MPI_File output_file;
	MPI_File_open(MPI_COMM_WORLD, output_filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
        MPI_File_set_view(output_file,
			  (MPI_Offset)metadata.offset,
			  MPI_UINT16_T,
			  writ_frame,
			  "native",
			  MPI_INFO_NULL);
	endian_swap(out_data, frame.writ.w * frame.writ.h);
	MPI_File_write(output_file,
		      out_data,
		      frame.writ.w * frame.writ.h,
		      MPI_UINT16_T,
		      MPI_STATUS_IGNORE);

	MPI_File_close(&output_file);
	MPI_Type_free(&writ_frame);	
	free(in_data);
	MPI_File_close(&input_file);
	MPI_Type_free(&read_frame);
	MPI_Type_free(&MPI_metadata);
	MPI_Finalize();
	return EXIT_SUCCESS;
}
