extern void composition(int dm[], double *A, double *B, double *C);
extern void composition_jacobian(int dm[],
                                 double *A, double * JA, double *B, double *JB,
                                 double *C, double *JC);
extern double samp(int dm[], double f[], double x, double y);
extern void expdef(int dm[], int k, double v[],
                   double t0[], double t1[], double J0[], double J1[]);
extern void unwrap(int dm[], double f[]);
extern int iteration_scratchsize(int dm[], int issym);
extern void iteration(int dm[], int k, double v[], double g[], double f[],
                      int rtype, double param[], double lmreg, int cycles, int nits, int issym, double ov[], double *buf);
extern void bracket(int dm[], double *A, double *B, double *C);

