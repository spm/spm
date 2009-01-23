/* $Id: optimizer3d.h 2645 2009-01-23 13:02:50Z john $ */
/* (c) John Ashburner (2007) */
extern void fmg3(int n0[], float *a0, float *b0, int rtype, double param[], int c, int nit,
                 float *u0, float *scratch);
extern void cgs3(int dm[], float A[], float b[], int rtype, double param[], double tol, int nit,
                 float x[], float r[], float p[], float Ap[]);
extern void resize(int na[], float *a,  int nc[], float *c, float *b);
extern float norm(int m, float a[]);
extern void LtLf_be(int dm[], float f[], double s[], float g[]);
extern void LtLf_me(int dm[], float f[], double s[], float g[]);
extern void LtLf_le(int dm[], float f[], double s[], float g[]);
extern int fmg3_scratchsize(int n0[]);

extern void fmg3_noa(int n0[], float *b0, int rtype, double param[], int c, int nit,
                 float *u0, float *scratch);
extern int fmg3_scratchsize_noa(int n0[]);

