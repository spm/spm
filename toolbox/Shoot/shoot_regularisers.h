/* $Id: shoot_regularisers.h 4583 2011-12-06 16:03:01Z john $ */
/* (c) John Ashburner (2011) */
extern void vel2mom(mwSize dm[], float f[], double s[], float g[]);
extern void relax(mwSize dm[], float a[], float b[], double s[], int nit, float u[]);
extern void Atimesp(mwSize dm[], float A[], double param[], float p[], float Ap[]);
extern double sumsq(mwSize dm[], float a[], float b[], double s[], float u[]);
