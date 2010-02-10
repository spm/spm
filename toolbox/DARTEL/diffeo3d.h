/* $Id: diffeo3d.h 3720 2010-02-10 18:26:58Z john $ */
/* (c) John Ashburner (2007) */
extern void composition(int dm[], float *A, float *B, float *C);
extern void composition_jacobian(int dm[],
                                 float *A, float * JA, float *B, float *JB,
                                 float *C, float *JC);
extern void composition_jacdet(int dm[],
                               float *A, float * JA, float *B, float *JB,
                               float *C, float *JC);
extern void smalldef(int dm[], double sc, float v[], float t[]);
extern void smalldef_jac1(int dm[], double sc, float v[], float t[], float J[]);
extern double samp(int dm[], float f[], double x, double y, double z);
extern void sampn(int dm[], float f[], int n, int mm, double x, double y, double z, double v[]);
extern void expdef(int dm[], int k, double sc, float v[],
                   float t0[], float t1[], float J0[], float J1[]);
extern void expdefdet(int dm[], int k, double sc, float v[],
                      float t0[], float t1[], float J0[], float J1[]);
extern void unwrap(int dm[], float f[]);
extern int iteration_scratchsize(int dm[], int issym, int k);
extern void iteration(int dm[], int k, float v[], float g[], float f[], float jd[],
                      int rtype, double param[], double lmreg0, int cycles, int its, int issym,
                      float ov[], double ll[], float *buf);
extern void bracket(int dm[], float *A, float *B, float *C);
extern void push(int dm[], int m, int n, float def[], float pf[], float po[], float so[]);
extern void pushc(int dm[], int m, int n, float def[], float pf[], float po[], float so[]);
extern void pushc_grads(int dm[], int m, float def[], float J[], float pf[], float po[]);
extern void determinant(int dm[], float J0[], float d[]);
extern void minmax_div(int dm[], float v0[], double mnmx[]);

