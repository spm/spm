/* $Id: shoot_regularisers.h 6799 2016-05-20 16:50:25Z john $ */
/* (c) John Ashburner (2011) */
extern double trapprox(mwSize dm[], float a[], double s[]);
extern void vel2mom(mwSize dm[], float f[], double s[], float g[]);
extern void relax(mwSize dm[], float a[], float b[], double s[], int nit, float u[]);
extern void Atimesp(mwSize dm[], float A[], double param[], float p[], float Ap[]);
extern double sumsq(mwSize dm[], float a[], float b[], double s[], float u[]);
extern void kernel(mwSize dm[], double s[], float f[]);

