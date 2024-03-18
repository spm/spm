/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#include "cuheader.h"
#define MAXVOL 20
#include "TVdenoise3d_dev.cu"

void TVdenoise3d(float *y, const float *x, const USIZE_t d[4], const float vox[3], const float lambdap[], const float lambdal[], const int nit)
{
    USIZE_t i, j, k;
    USIZE_t oi, oj, ok;
    int it;
    for(it=0; it<nit; it++)
        for(ok=0; ok<3; ok++)
            for(oj=0; oj<3; oj++)
                for(oi=0; oi<3; oi++)
                    for(k=1+ok; k<d[2]-1; k+=3)
                        for(j=1+oj; j<d[1]-1; j+=3)
                            for(i=1+oi; i<d[0]-1; i+=3)
                                TVdenoise3d_dev(i, j, k, y, x, d, vox, lambdap, lambdal);
}

void TVdenoise3d_fast(float *y, const float *x, const USIZE_t d[4], const float vox[3], const float lambdap[], const float lambdal[], const int nit)
{
    USIZE_t i, j, k;
    USIZE_t oi, oj, ok;
    int it;
    for(it=0; it<nit; it++)
        for(ok=0; ok<3; ok++)
            for(oj=0; oj<3; oj++)
                for(oi=0; oi<3; oi++)
                    for(k=1+ok; k<d[2]-1; k+=3)
                        for(j=1+oj; j<d[1]-1; j+=3)
                            for(i=1+oi; i<d[0]-1; i+=3)
                                TVdenoise3d_fast_dev(i, j, k, y, x, d, vox, lambdap, lambdal);
}
