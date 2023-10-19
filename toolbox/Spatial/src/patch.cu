/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#ifndef MAXD
#    define MAXD 5
#endif

__device__ void circ_b(SSIZE_t i, USIZE_t d, /*@out@*/USIZE_t *b, /*@out@*/int *s)
{
    *b = (USIZE_t)((i%(SSIZE_t)d) + d) % d;
    *s = 1;
}


__device__ void refl_b(SSIZE_t i, USIZE_t d, /*@out@*/USIZE_t *b, /*@out@*/int *s)
{
    USIZE_t d2 = d<<1;
    USIZE_t i1 = ((USIZE_t)(i%(SSIZE_t)d2) + d2) % d2;
    *b = (i1<d) ? i1 : d2-(i1+1);
    *s = 1;
}


__device__ void refl1_b(SSIZE_t i, USIZE_t d, /*@out@*/USIZE_t *b, /*@out@*/int *s)
{
    USIZE_t d2 = d<<1;
    USIZE_t i1 = ((USIZE_t)(i%(SSIZE_t)d2) + d2) % d2;
    if(i1<d)
    {
        *b = i1;
        *s = 1;
    }
    else
    {
        *b = d2-(i1+1);
        *s = -1;
    }
}


__device__ void bound(const int bnd, const SSIZE_t o, const USIZE_t dp, const USIZE_t d, /*@out@*/USIZE_t b[], /*@out@*/int s[])
{
    USIZE_t i;
    switch(bnd)
    {
    case 0:
        for(i=0; i<dp; i++)
            circ_b((SSIZE_t)i+o, d, b+i, s+i);
        break;
    case 1:
        for(i=0; i<dp; i++)
            refl_b((SSIZE_t)i+o, d, b+i, s+i);
        break;
    case 2:
        for(i=0; i<dp; i++)
            refl1_b((SSIZE_t)i+o, d, b+i, s+i);
        break;
    default:
        for(i=0; i<dp; i++)
            refl_b((SSIZE_t)i+o, d, b+i, s+i);
    }
}


__device__ void get_patch(const USIZE_t dp[], /*@out@*/float *fp, const int bnd[], SSIZE_t o[], const USIZE_t d0[], const float *f0)
{
    USIZE_t i, j, k;
    USIZE_t kb[MAXD], jb[MAXD], ib[MAXD];
    int    ks[MAXD], js[MAXD], is[MAXD];
    bound(bnd[0], o[0], dp[0], d0[0], ib, is);
    bound(bnd[1], o[1], dp[1], d0[1], jb, js);
    bound(bnd[2], o[2], dp[2], d0[2], kb, ks);

    for(k=0; k<dp[2]; k++)
    {
        USIZE_t tk  = d0[1]*kb[k];
        int    ksk = ks[k];
        for(j=0; j<dp[1]; j++)
        {
            USIZE_t tj = d0[0]*(jb[j] + tk);
            int    s  = ksk*js[j];
            for(i=0; i<dp[0]; i++, fp++)
                *fp = f0[ib[i] + tj]*(int)(s*is[i]);
        }
    }
}

__device__ float diffs(/*@out@*/float g[3], const int bnd[], SSIZE_t o[], const USIZE_t d[], const float *f0)
{
    USIZE_t b[3], t1, t2;
    int     s[3];
    float f;

    bound(bnd[0], o[0]-1, 3, d[0], b, s);
    t1   =        d[0]*(o[1] + d[1]*o[2]);
    g[0] = (f0[b[2]    + t1]*s[2] - f0[b[0] + t1]*s[0]);
    f    =  f0[o[0]    + t1];

    bound(bnd[1], o[1]-1, 3, d[1], b, s);
    t1   = o[0] + d[0]*(       d[1]*o[2]);
    t2   = d[0];
    g[1] = (f0[b[2]*t2 + t1]*s[2] - f0[b[0]*t2 + t1]*s[0]);

    bound(bnd[2], o[2]-1, 3, d[2], b, s);
    t1   = o[0] + d[0]* o[1];
    t2  *= d[1];
    g[2] = (f0[b[2]*t2 + t1]*s[2] - f0[b[0]*t2 + t1]*s[0]);

    return f;
}


__device__ void put_patch(const USIZE_t dp[], const float *fp, const int bnd[], const SSIZE_t o[], const USIZE_t d0[], float *f0)
{
    USIZE_t i, j, k;
    USIZE_t kb[MAXD], jb[MAXD], ib[MAXD];
    int    ks[MAXD], js[MAXD], is[MAXD];
    bound(bnd[0], o[0], dp[0], d0[0], ib, is);
    bound(bnd[1], o[1], dp[1], d0[1], jb, js);
    bound(bnd[2], o[2], dp[2], d0[2], kb, ks);

    for(k=0; k<dp[2]; k++)
    {
        USIZE_t tk  = d0[1]*kb[k];
        int    ksk = ks[k];
        for(j=0; j<dp[1]; j++)
        {
            USIZE_t tj = d0[0]*(jb[j] + tk);
            int    s  = ksk*js[j];
            for(i=0; i<dp[0]; i++, fp++)
            {
                float t = *fp*(int)(s*is[i]);
/*
#ifdef C
#               pragma omp atomic
#endif
*/
                atomicAdd(&(f0[ib[i] + tj]),t);
            }
        }
    }
}
