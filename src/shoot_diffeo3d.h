/* $Id: shoot_diffeo3d.h 7415 2018-09-10 18:24:16Z john $ */
/* (c) John Ashburner (2007) */
extern void composition(mwSize ma[], mwSize mm, float *A, float *B, float *C);
extern void composition_jacobian(mwSize ma[], mwSize mm,
                                 float *A, float * JA, float *B, float *JB,
                                 float *C, float *JC);
extern void composition_jacdet(mwSize dm[], mwSize mm,
                               float *A, float * JA, float *B, float *JB,
                               float *C, float *JC);
extern void smalldef(mwSize dm[], float sc, float v[], float t[]);
extern void smalldef_jac(mwSize dm[], float sc, float v0[], float t0[], float J0[]);
extern void smalldef_jac1(mwSize dm[], float sc, float v[], float t[], float J[]);
extern void sampn_vox(mwSize dm[], float f[], mwSize n, mwSize mm, double x, double y, double z, double v[]);
extern void unwrap(mwSize dm[], float f[]);
extern void bracket(mwSize dm[], float *A, float *B, float *C);
extern void push(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[], /*@out@*/ float F0[], /*@null@@out@*/ float S0[]);
extern void pushc(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[], /*@out@*/ float F0[], /*@null@@out@*/ float S0[]);
extern void pull(mwSize dm0[], mwSize m1, mwSize n, float Psi[],  float F0[], /*@out@*/ float F1[]);
extern void pullc(mwSize dm0[], mwSize m1, mwSize n, float Psi[],  float F0[], /*@out@*/ float F1[]);
extern void pushc_grads(mwSize dmo[], mwSize dmy[], float def[], float J[], float pf[], float po[]);
extern void determinant(mwSize dm[], float J0[], float d[]);
extern void minmax_div(mwSize dm[], float v0[], double mnmx[]);
extern void divergence(mwSize dm[], float v0[], float dv[]);
extern void def2det(mwSize dm[], float *Y, float *J, mwSignedIndex s);
extern void def2jac(mwSize dm[], float *Y, float *J, mwSignedIndex s);
extern void invdef(mwSize dim_y[3], float y[], mwSize dim_iy[3], float iy[], float M1[4][3], float M2[4][3]);
