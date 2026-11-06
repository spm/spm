/*
 * John Ashburner
 * Copyright (C) 2023 Department of Imaging Neuroscience, UCL
 */

#include "cuheader.h"
#define MAXVOL 20
#include "TVdenoise3d_dev.cu"

void TVdenoise3d(float *y, const float *x, const USIZE_t d[4], const float vox[3], const float lambdap[], const float lambdal[])
{
    USIZE_t i, j, k;
    USIZE_t oi, oj, ok;
    for(ok=0; ok<2; ok++)
        for(oj=0; oj<2; oj++)
            for(oi=0; oi<2; oi++)
                for(k=ok; k<d[2]; k+=2)
                    for(j=oj; j<d[1]; j+=2)
                        for(i=oi; i<d[0]; i+=2)
                            TVdenoise3d_dev(i, j, k, y, x, d, vox, lambdap, lambdal);
}

void TVdenoise3d_fast(float *y, const float *x, const USIZE_t d[4], const float vox[3], const float lambdap[], const float lambdal[])
{
    USIZE_t i, j, k;
    USIZE_t oi, oj, ok;
    for(ok=0; ok<2; ok++)
        for(oj=0; oj<2; oj++)
            for(oi=0; oi<2; oi++)
                for(k=ok; k<d[2]; k+=2)
                    for(j=oj; j<d[1]; j+=2)
                        for(i=oi; i<d[0]; i+=2)
                            TVdenoise3d_fast_dev(i, j, k, y, x, d, vox, lambdap, lambdal);
}

