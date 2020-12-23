#include <grid.h>
#include <math.h>

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
