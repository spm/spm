/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define MAXD 5
#include "cuheader.h"
#include "bspline.cu"
#include "patch.cu"

__device__ float dotp_patch(const USIZE_t dp[], const float fp[], const float wi[], const float wj[], const float wk[])
{
    USIZE_t i, j, k;
    float fo = 0.0f;
    for(k=0; k<dp[2]; k++)
    {
        float wkk = wk[k];
        USIZE_t ok = dp[1]*k;
        for(j=0; j<dp[1]; j++)
        {
            float wjj = wkk*wj[j];
            USIZE_t oj = dp[0]*(j+ok);
            for(i=0; i<dp[0]; i++)
            {
                fo += wi[i]*wjj*fp[oj+i];
            }
        }
    }
    return fo;
}


/*  Does the equivalent of:
    g[0] = dotp_patch(dp, fp, dwi,  wj,  wk);
    g[1] = dotp_patch(dp, fp,  wi, dwj,  wk);
    g[2] = dotp_patch(dp, fp,  wi,  wj, dwk); */
__device__ void grad_patch(const USIZE_t dp[], const float fp[],
                           const float   wi[], const float  wj[], const float  wk[],
                           const float  dwi[], const float dwj[], const float dwk[], /*@out@*/float g[])
{
    USIZE_t i, j, k;
    g[0] = g[1] = g[2] = 0.0f;

    for(k=0; k<dp[2]; k++)
    {
        float wkk = wk[k], wkkz = dwk[k];
        USIZE_t ok = dp[1]*k;
        for(j=0; j<dp[1]; j++)
        {
            float wjjx = wkk*wj[j], wjjy = wkk*dwj[j], wjjz = wkkz*wj[j];
            USIZE_t oj  = dp[0]*(j+ok);
            for(i=0; i<dp[0]; i++)
            {
                float tmp = fp[oj+i];
                g[0] += dwi[i]*wjjx*tmp;
                tmp  *= wi[i];
                g[1] += wjjy*tmp;
                g[2] += wjjz*tmp;
            }
        }
    }
}


__device__ float pull1(const USIZE_t d0[], const float f0[], const int bnd[], const USIZE_t dp[], const float x[])
{
    SSIZE_t o[3];
    float wi[MAXD], wj[MAXD], wk[MAXD],  fp[MAXD*MAXD*MAXD];
    o[0] = weights(dp[0], x[0], wi);
    o[1] = weights(dp[1], x[1], wj);
    o[2] = weights(dp[2], x[2], wk);
    get_patch(dp, fp, bnd, o, d0, f0);
    return dotp_patch(dp, fp, wi, wj, wk);
}


__device__ void pullg1(const USIZE_t d0[], const float f0[], const int bnd[], const USIZE_t dp[], const float x[], /*@out@*/float g[])
{
    SSIZE_t o[3];
    float  wi[MAXD],  wj[MAXD],  wk[MAXD];
    float dwi[MAXD], dwj[MAXD], dwk[MAXD];
    float fp[MAXD*MAXD*MAXD];
    o[0] = weights(dp[0], x[0],  wi);
    o[1] = weights(dp[1], x[1],  wj);
    o[2] = weights(dp[2], x[2],  wk);
    (void)dweights(dp[0], x[0], dwi);
    (void)dweights(dp[1], x[1], dwj);
    (void)dweights(dp[2], x[2], dwk);

    get_patch(dp, fp, bnd, o, d0, f0);
    grad_patch(dp, fp, wi, wj, wk, dwi, dwj, dwk, g);
}


__device__ void pullh1(const USIZE_t d0[], const float f0[], const int bnd[], const USIZE_t dp[], const float x[], /*@out@*/float h[])
{
    SSIZE_t o[3];
    float  wi[MAXD],  wj[MAXD],  wk[MAXD];
    float dwi[MAXD], dwj[MAXD], dwk[MAXD];
    float hwi[MAXD], hwj[MAXD], hwk[MAXD];
    float fp[MAXD*MAXD*MAXD];
    o[0] = weights(dp[0], x[0],  wi);
    o[1] = weights(dp[1], x[1],  wj);
    o[2] = weights(dp[2], x[2],  wk);
    (void)dweights(dp[0], x[0], dwi);
    (void)dweights(dp[1], x[1], dwj);
    (void)dweights(dp[2], x[2], dwk);
    (void)hweights(dp[0], x[0], hwi);
    (void)hweights(dp[1], x[1], hwj);
    (void)hweights(dp[2], x[2], hwk);

    get_patch(dp, fp, bnd, o, d0, f0);
    /* Slow. Could be speeded up, as in grad_patch */
    h[0] = dotp_patch(dp, fp, hwi,  wj,  wk);
    h[1] = dotp_patch(dp, fp,  wi, hwj,  wk);
    h[2] = dotp_patch(dp, fp,  wi,  wj, hwk);
    h[3] = dotp_patch(dp, fp, dwi, dwj,  wk);
    h[4] = dotp_patch(dp, fp, dwi,  wj, dwk);
    h[5] = dotp_patch(dp, fp,  wi, dwj, dwk);
}


__device__ void weight_patch(const USIZE_t dp[], /*@OUT@*/float fp[], const float wi[], const float wj[], const float wk[], const float fo)
{
    USIZE_t i, j, k;
    for(k=0; k<dp[2]; k++)
    {
        float wkk = wk[k];
        USIZE_t ok = dp[1]*k;
        for(j=0; j<dp[1]; j++)
        {
            float wjj = wkk*wj[j];
            USIZE_t oj = dp[0]*(j+ok);
            for(i=0; i<dp[0]; i++)
            {
                float wt = wi[i]*wjj;
                fp[oj+i] = wt*fo;
            }
        }
    }
}


__device__ void push1(const USIZE_t d0[], float f0[], /*@NULL@*/float c0[], const int bnd[], const USIZE_t dp[], const float x[], const float fo)
{
    SSIZE_t o[3];
    float wi[MAXD], wj[MAXD], wk[MAXD], fp[MAXD*MAXD*MAXD];
    o[0] = weights(dp[0], x[0], wi);
    o[1] = weights(dp[1], x[1], wj);
    o[2] = weights(dp[2], x[2], wk);
    weight_patch(dp, fp, wi, wj, wk, fo);
    put_patch(dp, fp, bnd, o, d0, f0);
    if(c0!=(void *)0)
    {
        weight_patch(dp, fp, wi, wj, wk, 1.0f);
        put_patch(dp, fp, bnd, o, d0, c0);
    }
}


__device__ void pushg1(const USIZE_t d0[], float f0[], const int bnd[], const USIZE_t dp[], const float x[], const float g[])
{
    SSIZE_t o[3];
    float  wi[MAXD],  wj[MAXD],  wk[MAXD];
    float dwi[MAXD], dwj[MAXD], dwk[MAXD];
    float fp[MAXD*MAXD*MAXD];
    o[0] = weights(dp[0], x[0],  wi);
    o[1] = weights(dp[1], x[1],  wj);
    o[2] = weights(dp[2], x[2],  wk);
    (void)dweights(dp[0], x[0], dwi);
    (void)dweights(dp[1], x[1], dwj);
    (void)dweights(dp[2], x[2], dwk);

    /* Slightly inefficient */
    weight_patch(dp, fp, dwi, wj, wk, g[0]);
    put_patch(dp, fp, bnd, o, d0, f0);

    weight_patch(dp, fp, wi, dwj, wk, g[1]);
    put_patch(dp, fp, bnd, o, d0, f0);

    weight_patch(dp, fp, wi, wj, dwk, g[2]);
    put_patch(dp, fp, bnd, o, d0, f0);
}


/* UNUSED */
__device__ void pushg1a(const USIZE_t d0[], float g0[], const int bnd[], const USIZE_t dp[], const float x[], const float f)
{
    SSIZE_t o[3], n0 = d0[0]*d0[1]*d0[2];
    float  wi[MAXD],  wj[MAXD],  wk[MAXD];
    float dwi[MAXD], dwj[MAXD], dwk[MAXD];
    float fp[MAXD*MAXD*MAXD];
    o[0] = weights(dp[0], x[0],  wi);
    o[1] = weights(dp[1], x[1],  wj);
    o[2] = weights(dp[2], x[2],  wk);
    (void)dweights(dp[0], x[0], dwi);
    (void)dweights(dp[1], x[1], dwj);
    (void)dweights(dp[2], x[2], dwk);

    /* Slightly inefficient */
    weight_patch(dp, fp, dwi, wj, wk, f);
    put_patch(dp, fp, bnd, o, d0, g0);

    weight_patch(dp, fp, wi, dwj, wk, f);
    put_patch(dp, fp, bnd, o, d0, g0+n0);

    weight_patch(dp, fp, wi, wj, dwk, f);
    put_patch(dp, fp, bnd, o, d0, g0+n0*2);
}


/* WORK IN PROGRESS */
__device__ void shootfun1_dev(const USIZE_t i, const USIZE_t j, const USIZE_t k, const USIZE_t d[], float u1[], const float u0[], const float v[], const int bnd[], const float s)
{
    SSIZE_t o[3];
    USIZE_t dp[3], n0 = d[0]*d[1]*d[2], oo = i+d[0]*(j+d[1]*k);
    float   x[3], uc[3], g[3],
	    J[3*3];    /* Gradients go into Jacobain matrix */

    float   s2 = -s*0.5f; /* account for two voxel spacing when computing gradients */

    o[0]     = i;
    o[1]     = j;
    o[2]     = k;

    x[0]     = diffs(g, bnd  , o, d, v     )*s + i;
    J[0+3*0] = s2*g[0] + 1.0f;
    J[0+3*1] = s2*g[1];
    J[0+3*2] = s2*g[2];

    x[1]     = diffs(g, bnd+3, o, d, v+n0  )*s + j;
    J[1+3*0] = s2*g[0];
    J[1+3*1] = s2*g[1] + 1.0f;
    J[1+3*2] = s2*g[2];

    x[2]     = diffs(g, bnd+6, o, d, v+n0*2)*s + k;
    J[2+3*0] = s2*g[0];
    J[2+3*1] = s2*g[1];
    J[2+3*2] = s2*g[2] + 1.0f;


    /* Get element of u0 at this point */
    uc[0] = u0[oo];
    uc[1] = u0[oo+n0];
    uc[2] = u0[oo+n0*2];

    /* Transform uc by Jacobain and push into u1 at position x (slow because of atomic additions) */
    dp[0] = dp[1] = dp[2] = 2; /* Trilinear */
    push1(d, u1     , (float *)NULL, bnd  , dp, x, J[0+3*0]*uc[0] + J[1+3*0]*uc[1] + J[2+3*0]*uc[2]);
    push1(d, u1+n0  , (float *)NULL, bnd+3, dp, x, J[0+3*1]*uc[0] + J[1+3*1]*uc[1] + J[2+3*1]*uc[2]);
    push1(d, u1+n0*2, (float *)NULL, bnd+6, dp, x, J[0+3*2]*uc[0] + J[1+3*2]*uc[1] + J[2+3*2]*uc[2]);
}


/* WORK IN PROGRESS */
__device__ void comp1_dev(const USIZE_t d[], float ff[], const float f0[], const int bnd[], const float x[])
{
    SSIZE_t o[3];
    USIZE_t n = d[0]*d[1]*d[2], dp[3];
    float wi[2], wj[2], wk[2],  fp[2*2*2];
    dp[0] = dp[1] = dp[2] = 2;
    o[0]  = weights(dp[0], x[0], wi);
    o[1]  = weights(dp[1], x[1], wj);
    o[2]  = weights(dp[2], x[2], wk);

    get_patch(dp, fp, bnd+3*0, o, d, f0    ); ff[0] = dotp_patch(dp, fp, wi, wj, wk);
    get_patch(dp, fp, bnd+3*1, o, d, f0+n  ); ff[1] = dotp_patch(dp, fp, wi, wj, wk);
    get_patch(dp, fp, bnd+3*2, o, d, f0+n*2); ff[2] = dotp_patch(dp, fp, wi, wj, wk);
}
