/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#include "cuheader.h"
#include "patch.cu"
#include "chol.cu"

__device__ float odconv(const float *lp, const float *vp, USIZE_t cp)
{
    const float *lpend;
    float vo   = 0.0f, vc = vp[cp];
    for(lpend=lp+cp; lp<lpend;) vo += *(lp++) * (*(vp++) - vc);
    vp++; /* Omit the central voxel */
    lp++;
    for(lpend=lp+cp; lp<lpend;) vo += *(lp++) * (*(vp++) - vc);
    return vo;
}


__device__ void relax1(USIZE_t i, USIZE_t j, USIZE_t k, float *v, const USIZE_t *d, const float *g, const float *h, const USIZE_t *dp, const float *l, const int *bnd)
{
    float vp[MAXD*MAXD*MAXD], b, a;
    SSIZE_t o[3];
    USIZE_t m, cp, hcode = d[4];
    o[0]  = (SSIZE_t)i-(SSIZE_t)dp[0]/2;
    o[1]  = (SSIZE_t)j-(SSIZE_t)dp[1]/2;
    o[2]  = (SSIZE_t)k-(SSIZE_t)dp[2]/2;

    cp    = (dp[0]*dp[1]*dp[2])/2;
    m     = i+d[0]*(j+d[1]*k);
    get_patch(dp, vp, bnd, o, d, v);
    b     = g[m] - odconv(l, vp, cp);
    a     = l[cp];
    if(hcode)
    {
        b -= h[m]*vp[cp];
        a += h[m];
    }
    v[m] += b/a;
}


#define MAXN 8

__device__ void relaxN(USIZE_t i, USIZE_t j, USIZE_t k, float *v, const USIZE_t *d, const float *g, const float *h, const USIZE_t *dp, const float *l, const int *bnd)
{
    float vp[MAXN*(MAXN+2)]; /* Also needs to big enough for prod(dp(1:3)) */
    SSIZE_t o[3];
    USIZE_t m, nd, np, cp, ii, d3 = d[3], hcode = d[4], lcode = dp[4];
    float *A = vp, *p = vp + MAXN*MAXN, *x = vp+MAXN*(MAXN+1); /* re-use memory */
    float b[MAXN];

    o[0] = i-(SSIZE_t)dp[0]/2;
    o[1] = j-(SSIZE_t)dp[1]/2;
    o[2] = k-(SSIZE_t)dp[2]/2;
    m    = i+d[0]*(j+d[1]*k);
    g   += m;
    h   += m;
    /* Original i, j & k no-longer needed */

    np   = dp[0]*dp[1]*dp[2];
    cp   = np/2;
    nd   = d[0]*d[1]*d[2];

    /* re-use i and j from here onwards */
    for(i=0; i<d3; i++)
        b[i] = g[nd*i];

    for(i=0, ii=d3; i<d3; i++)
    {
        USIZE_t j, ii0=ii;
        get_patch(dp, vp, bnd+i*3, o, d, v+i*nd);

        b[i] -= odconv(l+i*np, vp, cp);
        if(lcode==2)
        {
            for(j=i+1; j<d3; j++, ii++)
                b[j] -= 2.0f*odconv(l+ii*np, vp, cp);
        }
        if(hcode)
        {
            float vc = vp[cp];
            b[i] -= h[i*nd]*vc;

            if(hcode==2)
            {
                for(j=i+1, ii=ii0; j<d3; j++, ii++)
                    b[j] -= 2.0f*h[ii*nd]*vc;
            }
        }
    }

    /* Construct "diagonal" of L+H */
    l += cp;

    for(i=0; i<d3; i++)
        A[i+d3*i] = l[i*np]*1.000001f;
    if(lcode==2)
    {
        for(i=0, ii=d3; i<d3; i++)
        {
            USIZE_t j;
            for(j=i+1; j<d3; j++, ii++)
                A[i + j*d3] = A[j + i*d3] = l[ii*np];
        }
    }
    else
    {
        for(i=0, ii=d3; i<d3; i++)
        {
            USIZE_t j;
            for(j=i+1; j<d3; j++, ii++)
                A[i + j*d3] = A[j + i*d3] = 0.0f;
        }
    }

    if(hcode)
    {
        for(i=0; i<d3; i++)
            A[i+d3*i] += h[nd*i]*1.000001f;

        if(hcode==2)
        {
            for(i=0, ii=d3; i<d3; i++)
            {
                USIZE_t j;
                for(j=i+1; j<d3; j++, ii++)
                    A[i + j*d3] = A[j + i*d3] += h[ii*nd];
            }
        }
    }

    /* Compute x = A\b via Cholesky decomposition */
    choldcf(d3, A, p);
    chollsf(d3, A, p, b, x);

    v += m; /* shift pointer */
    for(i=0; i<d3; i++) v[nd*i] += x[i];
}

__device__ void convN(USIZE_t i, USIZE_t j, USIZE_t k, float *u, const float *v, const USIZE_t *d, const USIZE_t *dp, const float *l, const float *lc, const int *bnd)
{
    float vp[MAXD*MAXD*MAXD];
    SSIZE_t o[3];
    USIZE_t m, nd, np, cp, ii, d3 = d[3], lcode = dp[4];
    float b[MAXN];

    o[0] = i-(SSIZE_t)dp[0]/2;
    o[1] = j-(SSIZE_t)dp[1]/2;
    o[2] = k-(SSIZE_t)dp[2]/2;
    m    = i+d[0]*(j+d[1]*k);
    v   += m;
    np   = dp[0]*dp[1]*dp[2];
    cp   = np/2;
    nd   = d[0]*d[1]*d[2];

    /* re-use i and j from here onwards */
    for(i=0, ii=d3; i<d3; i++)
    {
        USIZE_t j;
        float vc = vp[cp];
        get_patch(dp, vp, bnd+i*3, o, d, v+i*nd);

        b[i] = odconv(l+i*np, vp, cp) + vc*lc[i];
        if(lcode==2)
        {
            for(j=i+1; j<d3; j++, ii++)
                b[j] += 2.0f*odconv(l+ii*np, vp, cp) + vc*lc[ii];
        }
    }

    u += m; /* shift pointer */
    for(i=0; i<d3; i++) u[nd*i] = b[i];
}
