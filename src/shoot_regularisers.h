/* $Id: shoot_regularisers.h 4875 2012-08-30 20:04:30Z john $ */
/* (c) John Ashburner (2011) */
extern void vel2mom(mwSize dm[], float f[], double s[], float g[]);
extern void relax(mwSize dm[], float a[], float b[], double s[], int nit, float u[]);
extern void Atimesp(mwSize dm[], float A[], double param[], float p[], float Ap[]);
extern double sumsq(mwSize dm[], float a[], float b[], double s[], float u[]);
extern void kernel(mwSize dm[], double s[], float f[]);
