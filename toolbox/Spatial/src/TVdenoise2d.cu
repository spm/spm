/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define CUDA
#include "cuheader.h"
#define SQUARE(x) (_t=(x), _t*_t)
#define SMALL 1e-5f

__device__ void TVdenoise2d_fast_dev(USIZE_t i, USIZE_t j, float y[], const float x[], const USIZE_t d[], const float lambda[])
{
    USIZE_t ij = i + j*d[0], n = d[0]*d[1], k;
    float _t, w22, w12 = 0.0f, w32 = 0.0f, w23 = 0.0f;
    float *yp, yb[12][4], *ybp;
    SSIZE_t d0 = (SSIZE_t)d[0];

    for(k=0, yp=y+ij; k<d[2]; k++, yp+=n)
    {
        float lambda_k = lambda[k];
        float y12, y13, y21, y22, y23, y31, y32;

        /* Doesn't handle edges, so need i>=1 & i<d[0]-1 & j>=1 & j<d[1]-1 */
        /* Use one of the four different arrangements for how the centre voxel (y22)
           could be influenced by its neighbours. Fast version just uses the
           first.
           o   o   o  :   o    y12   y13
                      :
               |   |  :         |     |
           o --* --o  :  y21 --y22 --y23
                      :
               |      :         |
           o --o   o  :  y31 --y32    o
        */
                           y12   = yp[-1]; y13   = yp[-1+d0];
        y21   = yp[  -d0]; y22   = *yp   ; y23   = yp[  +d0];
        y31   = yp[ 1-d0]; y32   = yp[ 1];

        ybp    = yb[k];
        ybp[0] = y12;
        ybp[1] = y32;
        ybp[2] = y21;
        ybp[3] = y23;

        w12  += lambda_k*(SQUARE(y22-y21) + SQUARE(y22-y12) + SMALL); /* links to 21 and 12 */
        w23  += lambda_k*(SQUARE(y22-y23) + SQUARE(y13-y23) + SMALL); /* link to 23 */
        w32  += lambda_k*(SQUARE(y22-y32) + SQUARE(y31-y32) + SMALL); /* link to 32 */
    }

    /* See https://francisbach.com/the-%ce%b7-trick-or-the-effectiveness-of-reweighted-least-squares */
    /* w.. = 1/eta.. */
    w12 = 1.0f/sqrt(w12);
    w32 = 1.0f/sqrt(w32);
    w23 = 1.0f/sqrt(w23);
    w22 = w12 + w12 + w32 + w23;
/*
    for(k=0, yp=y; k<d[2]; k++, yp+=n, x+=n)
        yp[ij] = (lambda[k]*((yp[ij-1] + yp[ij-0])*w12 + yp[ij+1]*w32 + yp[ij+d0]*w23) + x[ij])
                /(lambda[k]*w22+1.0f);
*/
    x += ij;
    for(k=0, yp=y+ij; k<d[2]; k++, x+=n, yp+=n)
    {
        ybp = yb[k];
        *yp = (lambda[k]*((ybp[0] + ybp[2])*w12 + ybp[1]*w32 + ybp[3]*w23) + *x)
             /(lambda[k]*w22+1.0f);
    }
}


__device__ void TVdenoise2d_dev(USIZE_t i, USIZE_t j, float y[], const float x[], const USIZE_t d[], const float lambda[])
{
    USIZE_t ij = i + j*d[0], n = d[0]*d[1], k;
    float _t, w22, w12, w32, w21, w23;
    float eta[12];
    float *yp, yb[12][4], *ybp;
    SSIZE_t d0 = (SSIZE_t)d[0];

    for(k=0; k<12; k++) eta[k] = 0.0f;

    for(k=0, yp=y+ij; k<d[2]; k++, yp+=n)
    {
        float lambda_k = lambda[k];
        float y11, y12, y13, y21, y22, y23, y31, y32, y33;
        float d12, d21, d32, d23;

        /* Doesn't handle edges, so need i>=1 & i<d[0]-1 & j>=1 & j<d[1]-1 */
        y11   = yp[-1-d0]; y12   = yp[-1]; y13   = yp[-1+d0];
        y21   = yp[  -d0]; y22   =*yp    ; y23   = yp[   d0];
        y31   = yp[ 1-d0]; y32   = yp[ 1]; y33   = yp[ 1+d0];

        /* Save for later */
        ybp    = yb[k];
        ybp[0] = y12;
        ybp[1] = y32;
        ybp[2] = y21;
        ybp[3] = y23;

        /* Four different arrangements for how the centre voxel (y22) 
           could be influenced by its neighbours.
         o   o   o  :  o   o-- o  :  o --o   o  :  o   o   o
                    :      |      :      |      :
             |   |  :             :             :  |   |
         o --* --o  :  o-- *-- o  :  o --* --o  :  o-- *-- o
                    :  |   |      :      |   |  :
             |      :             :             :      |
         o --o   o  :  o   o   o  :  o   o   o  :  o   o-- o
        */

        d12 = SQUARE(y22-y12);
        d21 = SQUARE(y22-y21);
        d23 = SQUARE(y22-y23);
        d32 = SQUARE(y22-y32);

        eta[0]  += lambda_k*(d21 + d12             + SMALL);
        eta[1]  += lambda_k*(d23 + SQUARE(y13-y23) + SMALL);
        eta[2]  += lambda_k*(d32 + SQUARE(y31-y32) + SMALL);

        eta[3]  += lambda_k*(d23 + d32             + SMALL);
        eta[4]  += lambda_k*(d21 + SQUARE(y31-y21) + SMALL);
        eta[5]  += lambda_k*(d12 + SQUARE(y13-y12) + SMALL);

        eta[6]  += lambda_k*(d32 + d21             + SMALL);
        eta[7]  += lambda_k*(d23 + SQUARE(y33-y23) + SMALL);
        eta[8]  += lambda_k*(d12 + SQUARE(y11-y12) + SMALL);

        eta[9]  += lambda_k*(d23 + d12             + SMALL);
        eta[10] += lambda_k*(d21 + SQUARE(y11-y21) + SMALL);
        eta[11] += lambda_k*(d32 + SQUARE(y33-y32) + SMALL);

    }

    /* See https://francisbach.com/the-%ce%b7-trick-or-the-effectiveness-of-reweighted-least-squares */
    for(k=0; k<12; k++) eta[k] = 1.0f/sqrt(eta[k]);

    w12 = (eta[0] + eta[5] + eta[8] + eta[ 9])*0.25f;
    w21 = (eta[0] + eta[4] + eta[6] + eta[10])*0.25f;
    w32 = (eta[2] + eta[3] + eta[6] + eta[11])*0.25f;
    w23 = (eta[1] + eta[3] + eta[7] + eta[ 9])*0.25f;

    w22 = w12 + w21 + w32 + w23;

    x  += ij;
    for(k=0, yp=y+ij; k<d[2]; k++, x+=n, yp+=n)
    {
        ybp = yb[k];
        *yp = (lambda[k]*(ybp[0]*w12 + ybp[1]*w32 + ybp[2]*w21 + ybp[3]*w23) + *x)
             /(lambda[k]*w22+1.0f);
    }
}



/* Use constant memory for lower access times */
__constant__ USIZE_t o[2];    /* Offsets (x & y) */
__constant__ USIZE_t d[3];    /* Image dimensions (x & y) and number of images*/
__constant__ float lambda[8]; /* Regularisation of each channel */

__global__ void TVdenoise2d(float *y, const float *x)
{
    USIZE_t i, j;

    /* Leaves edge voxels alone */
    j = (blockIdx.y*blockDim.y + threadIdx.y)*3 + 1 + o[1];
    if (j>=d[1]-1) return;

    i = (blockIdx.x*blockDim.x + threadIdx.x)*3 + 1 + o[0];
    if (i>=d[0]-1) return;

    TVdenoise2d_dev(i, j, y, x, d, lambda);
}

__global__ void TVdenoise2d_fast(float *y, const float *x)
{
    USIZE_t i, j;

    /* Leaves edge voxels alone */
    j = (blockIdx.y*blockDim.y + threadIdx.y)*3 + 1 + o[1];
    if (j>=d[1]-1) return;

    i = (blockIdx.x*blockDim.x + threadIdx.x)*3 + 1 + o[0];
    if (i>=d[0]-1) return;

    TVdenoise2d_fast_dev(i, j, y, x, d, lambda);
}
