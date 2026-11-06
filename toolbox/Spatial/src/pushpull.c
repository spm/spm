/*
 * John Ashburner
 * Copyright (C) 2023 Department of Imaging Neuroscience, UCL
 */

#include "cuheader.h"
#include "pushpull_dev.cu"
/* #define ISFINITE(x) (bool)(isfinite(x)) */
#define ISFINITE(x) isfinite(x)

#define VOXOK(x, d0) (ISFINITE(x[0]) && ISFINITE(x[1]) && ISFINITE(x[2]) && \
                      (ext!=0 || (x[0]>=-0.01f && x[0]<=(float)(d0[0])-0.99f && \
                                  x[1]>=-0.01f && x[1]<=(float)(d0[1])-0.99f && \
                                  x[2]>=-0.01f && x[2]<=(float)(d0[2])-0.99f)))



/* Pull voxels */
void pull(float *f1, const float *phi, const float *f0,
          const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext)
{
    USIZE_t i;
/*
#   pragma omp parallel for
*/
    for(i=0; i<n1; i++)
    {
        float x[3];
        x[0] = phi[i]       - 1.0f;
        x[1] = phi[i+n1]    - 1.0f;
        x[2] = phi[i+n1*2]  - 1.0f;

        if (VOXOK(x,d0))
            f1[i] = pull1(d0, f0, bnd, dp, x);
        else
            f1[i] = NAN;
    }
}


/* Pull gradients */
void pullg(float *g1, const float *phi, const float *f0,
          const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext)
{
    USIZE_t i;
/*
#   pragma omp parallel for
*/
    for(i=0; i<n1; i++)
    {
        float x[3];
        x[0] = phi[i]       - 1.0f;
        x[1] = phi[i+n1]    - 1.0f;
        x[2] = phi[i+n1*2]  - 1.0f;

        if (VOXOK(x,d0))
        {
            float g[3];
            pullg1(d0, f0, bnd, dp, x, g);
            g1[i       ] = g[0];
            g1[i + n1  ] = g[1];
            g1[i + n1*2] = g[2];
        }
        else
        {
            g1[i       ] = NAN;
            g1[i + n1  ] = NAN;
            g1[i + n1*2] = NAN;
        }
    }
}


/* Pull second derivatives */
void pullh(float *h1, const float *phi, const float *f0,
           const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext)
{
    USIZE_t i;
/*
#   pragma omp parallel for
*/
    for(i=0; i<n1; i++)
    {
        float x[3];
        x[0] = phi[i]        - 1.0f;
        x[1] = phi[i + n1]   - 1.0f;
        x[2] = phi[i + n1*2] - 1.0f;

        if (VOXOK(x,d0))
        {
            float h[6];
            pullh1(d0, f0, bnd, dp, x, h);
            h1[i       ] = h[0];
            h1[i + n1*4] = h[1];
            h1[i + n1*8] = h[2];
            h1[i + n1  ] = h1[i + n1*3] = h[3];
            h1[i + n1*2] = h1[i + n1*6] = h[4];
            h1[i + n1*5] = h1[i + n1*7] = h[5];
        }
        else
        {
            h1[i       ] = h1[i + n1  ] = h1[i + n1*2] = 
            h1[i + n1*3] = h1[i + n1*4] = h1[i + n1*5] =
            h1[i + n1*6] = h1[i + n1*7] = h1[i + n1*8] = NAN;
        }
    }
}


/* Push voxels */
void push(float *f0, const float *phi, const float *f1,
                  const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext)
{
    USIZE_t i;
    for(i=0; i<n1; i++)
    {
        float x[3], fo;
        fo   = f1[i];
        if(ISFINITE(fo))
        {
            x[0] = phi[i]       - 1.0f;
            x[1] = phi[i+n1]    - 1.0f;
            x[2] = phi[i+n1*2]  - 1.0f;

            if(VOXOK(x,d0)) push1(d0, f0, (float *)0, bnd, dp, x, fo);
        }
    }
}


/* Push gradients */
void pushg(float *g0, const float *phi, const float *g1,
           const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext)
{
    USIZE_t i;
    for(i=0; i<n1; i++)
    {
        float x[3], g[3];
        x[0] = phi[i]       - 1.0f;
        x[1] = phi[i+n1]    - 1.0f;
        x[2] = phi[i+n1*2]  - 1.0f;

        if (VOXOK(x,d0))
        {
            g[0] = g1[i];
            if(ISFINITE(g[0]))
            {
                g[1] = g1[i+n1];
                if(!ISFINITE(g[1]))
                {
                    g[2] = g1[i+n1*2];
                    if(ISFINITE(g[2])) pushg1(d0, g0, bnd, dp, x, g);
                }
            }
        }
    }
}


/* CURRENTLY UNUSED */
void pushg3(float *g0, const float *phi, const float *f1,
                   const USIZE_t *d0, const USIZE_t n1, const int *bnd, const USIZE_t *dp, const int ext)
{
    USIZE_t i;
    for(i=0; i<n1; i++)
    {
        float x[3], f;
        x[0] = phi[i]       - 1.0f;
        x[1] = phi[i+n1]    - 1.0f;
        x[2] = phi[i+n1*2]  - 1.0f;

        if (VOXOK(x,d0))
        {
            f = f1[i];
            if(ISFINITE(f)) pushg1a(d0, g0, bnd, dp, x, f);
        }
    }
}

/* WORK IN PROGRESS: Pull voxels */
void pull_affine(float *f1,         const float *f0,
                 const USIZE_t *d0, const USIZE_t n1,
                 const USIZE_t *d1, const float *Aff,
                 const int *bnd,    const USIZE_t *dp, const int ext)
{
    int i,j,k, ijk;
    /*
     * #   pragma omp parallel for
     */
    for(k=0, ijk=0; k<d1[2]; k++)
    {
        for(j=0; j<d1[1]; j++)
        {
            float xj[3];
            xj[0] = Aff[3]*j + Aff[6]*k + Aff[ 9];
            xj[1] = Aff[4]*j + Aff[7]*k + Aff[10];
            xj[2] = Aff[5]*j + Aff[8]*k + Aff[11];

            for(i=0; i<d1[0]; i++, ijk++)
            {
                float x[3];
                /* Assume Aff is adjusted for 0-offset data */
                x[0] = xj[0] + Aff[0]*i;
                x[1] = xj[1] + Aff[1]*i;
                x[2] = xj[2] + Aff[2]*i;

                if (VOXOK(x,d0))
                    f1[ijk] = pull1(d0, f0, bnd, dp, x);
                else
                    f1[ijk] = NAN;
            }
        }
    }
}


/* WORK IN PROGRESS */
void push_affine(float *f0,         const float *f1,
                 const USIZE_t *d0, const USIZE_t n1,
                 const USIZE_t *d1, const float *Aff,
                 const int *bnd,    const USIZE_t *dp, const int ext)
{
    int i,j,k, ijk;
    /* printf("[%d,%d,%d]\n", (int)(bnd[0]), (int)(bnd[1]), (int)(bnd[2])); */
    for(k=0, ijk=0; k<d1[2]; k++)
    {
        for(j=0; j<d1[1]; j++)
        {
            float xj[3];
            xj[0] = Aff[3]*j + Aff[6]*k + Aff[ 9];
            xj[1] = Aff[4]*j + Aff[7]*k + Aff[10];
            xj[2] = Aff[5]*j + Aff[8]*k + Aff[11];

            for(i=0; i<d1[0]; i++, ijk++)
            {
                float x[3], fo   = f1[ijk];
                if(ISFINITE(fo))
                {
                    /* Assume Aff is adjusted for 0-offset data */
                    x[0] = xj[0] + Aff[0]*i;
                    x[1] = xj[1] + Aff[1]*i;
                    x[2] = xj[2] + Aff[2]*i;

                    if(VOXOK(x,d0)) push1(d0, f0, (float *)0, bnd, dp, x, fo);
                }
            }
        }
    }
}

