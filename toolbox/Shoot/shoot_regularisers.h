/* $Id: shoot_regularisers.h 4573 2011-11-25 23:01:01Z john $ */
/* (c) John Ashburner (2011) */
extern void vel2mom_be(mwSize dm[], float f[], double s[], float g[]);
extern void vel2mom_me(mwSize dm[], float f[], double s[], float g[]);
extern void vel2mom_le(mwSize dm[], float f[], double s[], float g[]);

extern void relax_be(mwSize dm[], float a[], float b[], double s[], int nit, float u[]);
extern void relax_me(mwSize dm[], float a[], float b[], double s[], int nit, float u[]);
extern void relax_le(mwSize dm[], float a[], float b[], double s[], int nit, float u[]);

extern void Atimesp_be(mwSize dm[], float A[], double param[], float p[], float Ap[]);
extern void Atimesp_me(mwSize dm[], float A[], double param[], float p[], float Ap[]);
extern void Atimesp_le(mwSize dm[], float A[], double param[], float p[], float Ap[]);

