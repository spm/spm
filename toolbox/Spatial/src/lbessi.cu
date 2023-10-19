/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define CUDA
#include "cuheader.h"
#include "lbessi_dev.cu"

__device__ USIZE_t calculateGlobalIndex() {
    USIZE_t const globalBlockIndex = blockIdx.x + blockIdx.y * gridDim.x;
    USIZE_t const localThreadIdx   = threadIdx.x + blockDim.x * threadIdx.y;
    USIZE_t const threadsPerBlock  = blockDim.x*blockDim.y;
    return localThreadIdx + globalBlockIndex*threadsPerBlock;
}

__global__ void lbessi_element(float *out, const float nu, const float *z, const USIZE_t numel)
{
    USIZE_t const i = calculateGlobalIndex();
    if (i >= numel) return;
    out[i] = lbessif(nu,z[i]);
}
