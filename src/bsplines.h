/*
 * $Id: bsplines.h 7408 2018-08-24 14:54:57Z john $
 * John Ashburner
 */
#ifdef IMAGE_SINGLE
#define IMAGE_DTYPE float
#else
#define IMAGE_DTYPE double
#endif
 
void splinc_wrap(IMAGE_DTYPE c[], int m, double p[], int np);
void splinc_mirror(IMAGE_DTYPE c[], int m, double p[], int np);
int get_poles(int d, /*@out@*/ int *np, /*@out@*/ double p[]);
int mirror();
int wrap();
IMAGE_DTYPE sample3(IMAGE_DTYPE c[], int m0, int m1, int m2,
    IMAGE_DTYPE x0, IMAGE_DTYPE x1, IMAGE_DTYPE x2, int d[], int (*bnd[])());
IMAGE_DTYPE dsample3(IMAGE_DTYPE c[], int m0, int m1, int m2,
    IMAGE_DTYPE x0, IMAGE_DTYPE x1, IMAGE_DTYPE x2,
    int d[], IMAGE_DTYPE *pg0, IMAGE_DTYPE *pg1, IMAGE_DTYPE *pg2,
    int (*bnd[])());

