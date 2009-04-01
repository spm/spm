/* $Id: diffeo2d.h 3032 2009-04-01 14:14:18Z guillaume $ */
/* (c) John Ashburner (2007) */
extern void composition(int dm[], double *A, double *B, double *C);
extern void composition_jacobian(int dm[],
                                 double *A, double * JA, double *B, double *JB,
                                 double *C, double *JC);
extern void composition_detjac(int dm[],
                                 double *A, double * dA, double *B, double *dB,
                                 double *C, double *dC);
extern double samp(int dm[], double f[], double x, double y);
extern void expdef(int dm[], int k, double v[],
                   double t0[], double t1[], double J0[], double J1[]);
extern void expdefdet(int dm[], int k, double v[],
                   double t0[], double t1[], double J0[], double J1[]);
extern void unwrap(int dm[], double f[]);
extern int dartel_scratchsize(int dm[], int issym);
extern void dartel(int dm[], int k, double v[], double g[], double f[],
                   double dj[], int rtype, double param[], double lmreg,
                   int cycles, int nits, int issym,
                   double ov[], double ll[], double *buf);
extern void bracket(int dm[], double *A, double *B, double *C);

