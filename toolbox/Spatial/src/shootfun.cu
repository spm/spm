/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define CUDA
#include "cuheader.h"
#include "pushpull_dev.cu"
/* #include<math_constants.h> */

/* Use constant memory for lower access times */
__constant__ int     bnd[3*3]; /* boundary codes   */
__constant__ USIZE_t d[3];     /* image dimensions */
__constant__ float   s;        /* scaling of velocities */

__global__ void shootfun1_element(float *u1, const float *u0, const float *v0)
{
    USIZE_t i,j,k;

    i = (blockIdx.x*blockDim.x + threadIdx.x); if (i>=d[0]) return;
    j = (blockIdx.y*blockDim.y + threadIdx.y); if (j>=d[1]) return;
    k = (blockIdx.z*blockDim.z + threadIdx.z); if (k>=d[2]) return;

    shootfun1_dev(i, j, k, d, u1, u0, v0, bnd, s);
}

__global__ void shootfun2_element(float *psi1, const float *psi0, const float *v)
{
    USIZE_t i,j,k,o,n;
    float x[3], ff[3];

    i = (blockIdx.x*blockDim.x + threadIdx.x); if (i>=d[0]) return;
    j = (blockIdx.y*blockDim.y + threadIdx.y); if (j>=d[1]) return;
    k = (blockIdx.z*blockDim.z + threadIdx.z); if (k>=d[2]) return;
    o = i+d[0]*(j+d[1]*k);
    n = d[0]*d[1]*d[2];

    x[0] = (float)i-v[o+n*0]*s;
    x[1] = (float)j-v[o+n*1]*s;
    x[2] = (float)k-v[o+n*2]*s;
    comp1_dev(d, ff, psi0, bnd, x);
    psi1[o+n*0] = ff[0]-v[o+n*0]*s;
    psi1[o+n*1] = ff[1]-v[o+n*1]*s;
    psi1[o+n*2] = ff[2]-v[o+n*2]*s;
}
