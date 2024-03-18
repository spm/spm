/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#include "cuheader.h"
#include "pushpull_dev.cu"

void shootfun1(float *u1, const float *u0, const float *v0, const int bnd[], const USIZE_t d[], const float s)
{
    USIZE_t i,j,k;
    for(k=0; k<d[2]; k++)
        for(j=0; j<d[1]; j++)
            for(i=0; i<d[0]; i++)
                shootfun1_dev(i, j, k, d, u1, u0, v0, bnd, s);
}

void shootfun2(float *psi1, const float *psi0, const float *v, const int bnd[], const USIZE_t d[], const float s)
{
    USIZE_t i,j,k;

    for(k=0; k<d[2]; k++)
        for(j=0; j<d[1]; j++)
            for(i=0; i<d[0]; i++)
            {
                USIZE_t o = i+d[0]*(j+d[1]*k),
                        n = d[0]*d[1]*d[2];
                float x[3], ff[3];

                x[0] = (float)i-v[o+n*0]*s;
                x[1] = (float)j-v[o+n*1]*s;
                x[2] = (float)k-v[o+n*2]*s;
                comp1_dev(d, ff, psi0, bnd, x);
                psi1[o+n*0] = ff[0]-v[o+n*0]*s;
                psi1[o+n*1] = ff[1]-v[o+n*1]*s;
                psi1[o+n*2] = ff[2]-v[o+n*2]*s;
            }
}
