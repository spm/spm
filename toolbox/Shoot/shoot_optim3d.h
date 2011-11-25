/* $Id: shoot_optim3d.h 4573 2011-11-25 23:01:01Z john $ */
/* (c) John Ashburner (2011) */
extern void fmg3(mwSize n0[], float *a0, float *b0, int rtype, double param[], int c, int nit,
                 float *u0, float *scratch);
extern void cgs3(mwSize dm[], float A[], float b[], int rtype, double param[], double tol, int nit,
                 float x[], float r[], float p[], float Ap[]);
extern void resize(mwSize na[], float *a, mwSize nc[], float *c, float *b);
extern void restrict_vol(mwSize na[], float *a, mwSize nc[], float *c, float *b);
extern float norm(mwSize m, float a[]);
extern void vel2mom_be(mwSize dm[], float f[], double s[], float g[]);
extern void vel2mom_me(mwSize dm[], float f[], double s[], float g[]);
extern void vel2mom_le(mwSize dm[], float f[], double s[], float g[]);
extern mwSize fmg3_scratchsize(mwSize n0[], int use_hessian);

