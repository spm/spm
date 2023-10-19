/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define USIZE_t unsigned long long /* Done because of problems in Windows */

extern void      TVdenoise3d(float *y, const float *x, const USIZE_t d[4], const float vox[3], const float lambdap[], const float lambdal[], const int nit);
extern void TVdenoise3d_fast(float *y, const float *x, const USIZE_t d[4], const float vox[3], const float lambdap[], const float lambdal[], const int nit);
