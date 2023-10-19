/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define CUDA
#include "cuheader.h"
#include "operator_dev.cu"
/*
    MAXD - maximum filter size
    MAXN - maximum number of gradient fields
*/

/*
    Use constant memory for faster access times
    Limited by available constant memory on device:
        CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY
*/
__constant__ float   l[MAXD*MAXD*MAXD*MAXN*(MAXN+1)/2];  /* Filters */
__constant__ float   lc[MAXN*(MAXN+1)/2]; /* sum over filter elements */
__constant__ int     bnd[3*MAXN];         /* Boundary conditions */
__constant__ USIZE_t dp[5];               /* filter dimensions */
__constant__ USIZE_t  o[3];               /* offsets into volume */
__constant__ USIZE_t  d[5];               /* image data dimensions */
__constant__ USIZE_t  n[3];               /* number of elements */


__global__ void relax_element(float *v, const float *g, const float *h)
{
    USIZE_t ijk, i, j, k;
    ijk = threadIdx.x + blockDim.x*blockIdx.x;
    i   = (ijk % n[0])*dp[0] + o[0]; if(i>=d[0]) return;
    ijk =  ijk / n[0];
    j   = (ijk % n[1])*dp[1] + o[1]; if(j>=d[1]) return;
    k   = (ijk / n[1])*dp[2] + o[2]; if(k>=d[2]) return;

    relaxN(i, j, k, v, d, g, h, dp, l, bnd);
}


__global__ void conv_element(float *u, const float *v)
{
    USIZE_t i, j, k;
    i = threadIdx.x + blockDim.x*blockIdx.x; if(i>=d[0]) return;
    j = threadIdx.y + blockDim.y*blockIdx.y; if(j>=d[1]) return;
    k = threadIdx.z + blockDim.z*blockIdx.z; if(k>=d[2]) return;

    convN(i, j, k, u, v, d, dp, l, lc, bnd);
}
