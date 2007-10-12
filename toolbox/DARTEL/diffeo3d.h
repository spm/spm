/* (c) John Ashburner (2007) */
extern void composition(int dm[], float *A, float *B, float *C);
extern void composition_jacobian(int dm[],
                                 float *A, float * JA, float *B, float *JB,
                                 float *C, float *JC);
extern void composition_jacdet(int dm[],
                               float *A, float * JA, float *B, float *JB,
                               float *C, float *JC);

extern double samp(int dm[], float f[], double x, double y, double z);
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

