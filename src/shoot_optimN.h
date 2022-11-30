/*
 * John Ashburner
 * Copyright (C) 2011-2022 Wellcome Centre for Human Neuroimaging
 */

extern void fmg(mwSize n0[], float *a0, float *b0, double param[], double scal[], int c, int nit, float *u0, float *scratch);
extern void solve(mwSize dm[], float a[], float b[], double s[], double scal[], float u[]);
extern void LtLf(mwSize dm[], float f[], double s[], double scal[], float g[]);
extern void LtWLf(mwSize dm[], float f[], float h[], double s[], double scal[], float g[]);
extern void Atimesp1(mwSize dm[], float A[], float p[], float Ap[]);
extern double diaginv(mwSize dm[], float a[], float b[], double s[], double scal[], float u[]);
extern mwSize fmg_scratchsize(mwSize n0[]);
