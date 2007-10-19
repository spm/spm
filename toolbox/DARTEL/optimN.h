/* $Id: optimN.h 964 2007-10-19 16:35:34Z john $ */
/* (c) John Ashburner (2007) */
extern void fmg(int n0[], float *a0, float *b0, int rtype, double param[], int c, int nit,
                 float *u0, float *scratch);
extern void resize(int na[], float *a, int nc[], float *c, float *b);
extern float norm(int m, float a[]);
extern void LtLf_be(int dm[], float f[], double s[], float g[]);
extern void LtLf_me(int dm[], float f[], double s[], float g[]);
extern int fmg_scratchsize(int n0[]);

