/* $Id: shoot_optimN.h 7464 2018-10-31 16:57:27Z john $ */
/* (c) John Ashburner (2007) */
extern void fmg(mwSize n0[], float *a0, float *b0, double param[], double scal[], int c, int nit, float *u0, float *scratch);
extern void solve(mwSize dm[], float a[], float b[], double s[], double scal[], float u[]);
extern void LtLf(mwSize dm[], float f[], double s[], double scal[], float g[]);
extern int fmg_scratchsize(mwSize n0[]);

