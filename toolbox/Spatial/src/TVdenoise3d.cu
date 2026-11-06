/*
 * John Ashburner
 * Copyright (C) 2023 Department of Imaging Neuroscience, UCL
 */

#define CUDA
#include "cuheader.h"
#define MAXVOL 20
#include "TVdenoise3d_dev.cu"

/* Use constant memory for lower access times */
__constant__ USIZE_t o[3];    /* Offsets (x & y) */
__constant__ float vox[3];    /* Voxel sizes */
__constant__ USIZE_t d[4];    /* Image dimensions (x & y) and number of images*/
__constant__ float lambdap[MAXVOL]; /* Regularisation of each channel */
__constant__ float lambdal[MAXVOL]; /* Reciprocals of variances of each channel */

__global__ void TVdenoise3d(float *y, const float *x)
{
    USIZE_t i, j, k;

    /* Leaves edge voxels alone */
    k = (blockIdx.z*blockDim.z + threadIdx.z)*2 + o[2];
    if (k>=d[2]) return;

    j = (blockIdx.y*blockDim.y + threadIdx.y)*2 + o[1];
    if (j>=d[1]) return;

    i = (blockIdx.x*blockDim.x + threadIdx.x)*2 + o[0];
    if (i>=d[0]) return;

    TVdenoise3d_dev(i, j, k, y, x, d, vox, lambdap, lambdal);
}

__global__ void TVdenoise3d_fast(float *y, const float *x)
{
    USIZE_t i, j, k;

    /* Leaves edge voxels alone */
    k = (blockIdx.z*blockDim.z + threadIdx.z)*2 + o[2];
    if (k>=d[2]) return;

    j = (blockIdx.y*blockDim.y + threadIdx.y)*2 + o[1];
    if (j>=d[1]) return;

    i = (blockIdx.x*blockDim.x + threadIdx.x)*2 + o[0];
    if (i>=d[0]) return;

    TVdenoise3d_fast_dev(i, j, k, y, x, d, vox, lambdap, lambdal);
}

