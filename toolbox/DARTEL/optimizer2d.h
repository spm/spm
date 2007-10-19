/* $Id: optimizer2d.h 964 2007-10-19 16:35:34Z john $ */
/* (c) John Ashburner (2007) */
extern void fmg2(int n0[], double *a0, double *b0, int rtype, double param[], int c, int nit,
                 double *u0, double *scratch);
extern void cgs2(int dm[], double A[], double b[], int rtype, double param[], double tol, int nit,
                 double x[], double r[], double p[], double Ap[]);
extern void resize(int na[], double *a, int nc[], double *c, double *b);
extern double norm(int m, double a[]);
extern void LtLf_be(int dm[], double f[], double s[], double g[]);
extern void LtLf_me(int dm[], double f[], double s[], double g[]);
extern void LtLf_le(int dm[], double f[], double s[], double g[]);
extern int fmg2_scratchsize(int n0[]);

