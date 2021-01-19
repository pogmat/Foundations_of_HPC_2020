#ifndef _GRID_H
#define _GRID_H

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

static inline int min(const int a, const int b) {return((a < b) ? a : b);}
static inline int max(const int a, const int b) {return((a > b) ? a : b);}

void grid_dimension(const int, const int,
		    const int,
		    int* const, int* const);
void get_frame(const int, const int,
	       const int, const int,
	       const int, const int,
	       const int,
	       thimg_t* const);

#endif
