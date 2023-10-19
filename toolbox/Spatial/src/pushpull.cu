/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define CUDA
#include "cuheader.h"
#include "pushpull_dev.cu"
#include<math_constants.h>
#define ISFINITE(x) isfinite(x)

__device__ USIZE_t calculateGlobalIndex() {
    USIZE_t const globalBlockIndex = blockIdx.x + blockIdx.y * gridDim.x;
    USIZE_t const localThreadIdx   = threadIdx.x + blockDim.x * threadIdx.y;
    USIZE_t const threadsPerBlock  = blockDim.x*blockDim.y;
    return localThreadIdx + globalBlockIndex*threadsPerBlock;
}


/* Use constant memory for lower access times */
__constant__ int     bnd[3];  /* boundary codes   */
__constant__ USIZE_t  dp[3];  /* patch dimensions */
__constant__ USIZE_t  d0[3];  /* image dimensions */
__constant__ USIZE_t  n1;     /* Number of voxels */
__constant__ int     ext;     /* Extrapolate flag */

/* These are used by affine push/pull */
__constant__ USIZE_t  d1[3];  /* image dimensions */
__constant__ float   Aff[12]; /* Part of affine transform */


#define VOXOK(x, d0) (ISFINITE(x[0]) && ISFINITE(x[1]) && ISFINITE(x[2]) && \
                      (ext!=0 || (x[0]>=-0.01f && x[0]<=(float)(d0[0])-0.99f && \
                                  x[1]>=-0.01f && x[1]<=(float)(d0[1])-0.99f && \
                                  x[2]>=-0.01f && x[2]<=(float)(d0[2])-0.99f)))


__global__ void pull_element(float *f1, const float *phi, const float *f0)
{
    USIZE_t const globalThreadIdx = calculateGlobalIndex();
    float x[3];

    if(globalThreadIdx >= n1)
        return;

    x[0] = phi[globalThreadIdx]       - 1.0f;
    x[1] = phi[globalThreadIdx+n1]    - 1.0f;
    x[2] = phi[globalThreadIdx+n1*2]  - 1.0f;

    if(VOXOK(x,d0))
        f1[globalThreadIdx] = pull1(d0, f0, bnd, dp, x);
    else
        f1[globalThreadIdx] = CUDART_NAN_F;
}


__global__ void pullg_element(float *g1, const float *phi, const float *f0)
{
    USIZE_t const globalThreadIdx = calculateGlobalIndex();
    float x[3];

    if(globalThreadIdx >= n1) return;

    x[0] = phi[globalThreadIdx]       - 1.0f;
    x[1] = phi[globalThreadIdx+n1]    - 1.0f;
    x[2] = phi[globalThreadIdx+n1*2]  - 1.0f;

    if(VOXOK(x,d0))
    {
        float g[3];
        pullg1(d0, f0, bnd, dp, x, g);
        g1[globalThreadIdx       ] = g[0];
        g1[globalThreadIdx + n1  ] = g[1];
        g1[globalThreadIdx + n1*2] = g[2];
    }
    else
    {
        g1[globalThreadIdx       ] = CUDART_NAN_F;
        g1[globalThreadIdx + n1  ] = CUDART_NAN_F;
        g1[globalThreadIdx + n1*2] = CUDART_NAN_F;
    }
}


__global__ void pullh_element(float *h1, const float *phi, const float *f0)
{
    USIZE_t const globalThreadIdx = calculateGlobalIndex();
    float x[3];

    if(globalThreadIdx >= n1) return;

    x[0] = phi[globalThreadIdx]       - 1.0f;
    x[1] = phi[globalThreadIdx+n1]    - 1.0f;
    x[2] = phi[globalThreadIdx+n1*2]  - 1.0f;

    if(VOXOK(x,d0))
       {
           float h[6];
           pullh1(d0, f0, bnd, dp, x, h);
           h1[globalThreadIdx       ] = h[0];
           h1[globalThreadIdx + n1*4] = h[1];
           h1[globalThreadIdx + n1*8] = h[2];
           h1[globalThreadIdx + n1  ] = h1[globalThreadIdx + n1*3] = h[3];
           h1[globalThreadIdx + n1*2] = h1[globalThreadIdx + n1*6] = h[4];
           h1[globalThreadIdx + n1*5] = h1[globalThreadIdx + n1*7] = h[5];
       }
       else
       {
           h1[globalThreadIdx       ] = h1[globalThreadIdx + n1  ] = h1[globalThreadIdx + n1*2] = 
           h1[globalThreadIdx + n1*3] = h1[globalThreadIdx + n1*4] = h1[globalThreadIdx + n1*5] = 
           h1[globalThreadIdx + n1*6] = h1[globalThreadIdx + n1*7] = h1[globalThreadIdx + n1*8] = CUDART_NAN_F;
       }
}


__global__ void push_element(float *f0, const float *phi, const float *f1)
{
    USIZE_t const globalThreadIdx = calculateGlobalIndex();
    float x[3];
    float fo;

    if(globalThreadIdx >= n1) return;
    fo   = f1[globalThreadIdx];               if(!isfinite(fo))   return;
    x[0] = phi[globalThreadIdx]       - 1.0f; if(!isfinite(x[0])) return;
    x[1] = phi[globalThreadIdx+n1]    - 1.0f; if(!isfinite(x[1])) return;
    x[2] = phi[globalThreadIdx+n1*2]  - 1.0f; if(!isfinite(x[2])) return;

    if(ext || (x[0]>=-0.01f && x[0]<=d0[0]-0.99f &&
                x[1]>=-0.01f && x[1]<=d0[1]-0.99f &&
                x[2]>=-0.01f && x[2]<=d0[2]-0.99f))
        push1(d0, f0, (float *)0, bnd, dp, x, fo);
}


__global__ void pushc_element(float *f0, float *c0, const float *phi, const float *f1)
{
    USIZE_t const globalThreadIdx = calculateGlobalIndex();
    float x[3];
    float fo;

    if(globalThreadIdx >= n1) return;
    fo   = f1[globalThreadIdx];               if(!isfinite(fo))   return;
    x[0] = phi[globalThreadIdx]       - 1.0f; if(!isfinite(x[0])) return;
    x[1] = phi[globalThreadIdx+n1]    - 1.0f; if(!isfinite(x[1])) return;
    x[2] = phi[globalThreadIdx+n1*2]  - 1.0f; if(!isfinite(x[2])) return;

    if(ext || (x[0]>=-0.01f && x[0]<=d0[0]-0.99f &&
               x[1]>=-0.01f && x[1]<=d0[1]-0.99f &&
               x[2]>=-0.01f && x[2]<=d0[2]-0.99f))
        push1(d0, f0, c0, bnd, dp, x, fo);
}


__global__ void pushg_element(float *f0, const float *phi, const float *g1)
{
    USIZE_t const globalThreadIdx = calculateGlobalIndex();
    float x[3];
    float g[3];

    if(globalThreadIdx >= n1) return;
    x[0] = phi[globalThreadIdx]       - 1.0f;
    x[1] = phi[globalThreadIdx+n1]    - 1.0f;
    x[2] = phi[globalThreadIdx+n1*2]  - 1.0f;

    if(VOXOK(x,d0))
    {
        g[0] = g1[globalThreadIdx];
        g[1] = g1[globalThreadIdx + n1];
        g[2] = g1[globalThreadIdx + n1*2];
        pushg1(d0, f0, bnd, dp, x, g);
    }
}


/* UNUSED */
__global__ void pushg3_element(float *g0, const float *phi, const float *f1)
{
    USIZE_t const globalThreadIdx = calculateGlobalIndex();
    float x[3];
    float f;

    if(globalThreadIdx >= n1) return;
    x[0] = phi[globalThreadIdx]       - 1.0f; if(!isfinite(x[0])) return;
    x[1] = phi[globalThreadIdx+n1]    - 1.0f; if(!isfinite(x[1])) return;
    x[2] = phi[globalThreadIdx+n1*2]  - 1.0f; if(!isfinite(x[2])) return;

    if(ext || (x[0]>=-0.01f && x[0]<=d0[0]-0.99f &&
                x[1]>=-0.01f && x[1]<=d0[1]-0.99f &&
                x[2]>=-0.01f && x[2]<=d0[2]-0.99f))
    {
        f = f1[globalThreadIdx];  if(!ISFINITE(f)) return;
        pushg1a(d0, g0, bnd, dp, x, f);
    }
}



/* WORK IN PROGRESS */
__global__ void affine_pull_element(float *f1, const float *f0)
{
    USIZE_t globalThreadIdx;
    int i,j,k; 
    float x[3];
    
    i = (int)(blockIdx.x*blockDim.x + threadIdx.x); if (i>=d1[0]) return;
    j = (int)(blockIdx.y*blockDim.y + threadIdx.y); if (j>=d1[1]) return;
    k = (int)(blockIdx.z*blockDim.z + threadIdx.z); if (k>=d1[2]) return;

    globalThreadIdx = (USIZE_t)i + d1[0]*((USIZE_t)j + d1[1]*(USIZE_t)k);

    /* Assume Aff is adjusted for 0-offset data */
    x[0] = Aff[0]*i + Aff[3]*j + Aff[6]*k + Aff[ 9];
    x[0] = Aff[1]*i + Aff[4]*j + Aff[7]*k + Aff[10];
    x[0] = Aff[2]*i + Aff[5]*j + Aff[8]*k + Aff[11];

    if(VOXOK(x,d0))
        f1[globalThreadIdx] = pull1(d0, f0, bnd, dp, x);
    else
        f1[globalThreadIdx] = CUDART_NAN_F;
}


/* WORK IN PROGRESS */
__global__ void affine_push_element(float *f0, const float *f1)
{
    USIZE_t globalThreadIdx;
    int i,j,k;
    float x[3];
    float fo;

    i = (int)(blockIdx.x*blockDim.x + threadIdx.x); if (i>=d1[0]) return;
    j = (int)(blockIdx.y*blockDim.y + threadIdx.y); if (j>=d1[1]) return;
    k = (int)(blockIdx.z*blockDim.z + threadIdx.z); if (k>=d1[2]) return;

    globalThreadIdx = (USIZE_t)i + d1[0]*((USIZE_t)j + d1[1]*(USIZE_t)k);
    fo   = f1[globalThreadIdx]; if(!isfinite(fo)) return;

    /* Assume Aff is adjusted for 0-offset data */
    x[0] = Aff[0]*i + Aff[3]*j + Aff[6]*k + Aff[ 9];
    x[0] = Aff[1]*i + Aff[4]*j + Aff[7]*k + Aff[10];
    x[0] = Aff[2]*i + Aff[5]*j + Aff[8]*k + Aff[11];

    if(ISFINITE(x[0]) && ISFINITE(x[1]) && ISFINITE(x[2]) &&
       (ext || (x[0]>=-0.01f && x[0]<=d0[0]-0.99f &&
                x[1]>=-0.01f && x[1]<=d0[1]-0.99f &&
                x[2]>=-0.01f && x[2]<=d0[2]-0.99f)))
        push1(d0, f0, (float *)0, bnd, dp, x, fo);
}
