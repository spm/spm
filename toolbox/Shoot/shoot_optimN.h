/* $Id: optimN.h 2645 2009-01-23 13:02:50Z john $ */
/* (c) John Ashburner (2007) */
extern void fmg(mwSize n0[], float *a0, float *b0, double param[], double scal[], int c, int nit,
                 float *u0, float *scratch);
extern void LtLf(mwSize dm[], float f[], double s[], double scal[], float g[]);
extern int fmg_scratchsize(mwSize n0[]);

