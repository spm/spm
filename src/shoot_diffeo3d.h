/* $Id: shoot_diffeo3d.h 4875 2012-08-30 20:04:30Z john $ */
/* (c) John Ashburner (2007) */
extern void composition(mwSize ma[], mwSize mm, float *A, float *B, float *C);
extern void composition_jacobian(mwSize ma[], mwSize mm,
                                 float *A, float * JA, float *B, float *JB,
                                 float *C, float *JC);
extern void composition_jacdet(mwSize dm[], mwSize mm,
                               float *A, float * JA, float *B, float *JB,
                               float *C, float *JC);
extern void smalldef(mwSize dm[], double sc, float v[], float t[]);
extern void smalldef_jac(mwSize dm[], double sc, float v0[], float t0[], float J0[]);
extern void smalldef_jac1(mwSize dm[], double sc, float v[], float t[], float J[]);
extern double samp(mwSize dm[], float f[], double x, double y, double z);
extern void sampn(mwSize dm[], float f[], mwSize n, mwSize mm, double x, double y, double z, double v[]);
extern void unwrap(mwSize dm[], float f[]);
extern void bracket(mwSize dm[], float *A, float *B, float *C);
extern void push(mwSize dm[], mwSize m, mwSize n, float def[], float pf[], float po[], float so[]);
extern void pushc(mwSize dm[], mwSize m, mwSize n, float def[], float pf[], float po[], float so[]);
extern void pushc_grads(mwSize dmo[], mwSize dmy[], float def[], float J[], float pf[], float po[]);
extern void determinant(mwSize dm[], float J0[], float d[]);
extern void minmax_div(mwSize dm[], float v0[], double mnmx[]);
extern void divergence(mwSize dm[], float v0[], float dv[]);
extern void def2det(mwSize dm[], float *Y, float *J, mwSignedIndex s);
extern void def2jac(mwSize dm[], float *Y, float *J, mwSignedIndex s);
extern void invdef(mwSize dim_y[3], float y[], mwSize dim_iy[3], float iy[], double M1[4][3], double M2[4][3]);
