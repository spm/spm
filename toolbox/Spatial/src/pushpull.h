/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define USIZE_t unsigned long long /* Done because of problems in Windows */

extern void  push(float *f0, const float *phi, const float *f1,
                  const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext);
extern void pushg(float *f0, const float *phi, const float *f1,
                  const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext);
extern void  pull(float *f1, const float *phi, const float *f0,
                  const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext);
extern void pullg(float *g1, const float *phi, const float *f0,
                  const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext);
/*
extern void shootfun1(float *u1, const float *u0, const float *v0, const int bnd[], const USIZE_t d[], const float s);

extern void shootfun2(float *psi1, const float *psi0, const float *v, const int bnd[], const USIZE_t d[], const float s);
*/
