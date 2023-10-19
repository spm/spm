/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#include "cuheader.h"

/* Compute log(bessi(nu, z))
   - requires z to be reasonably small to work efficiently
*/
__device__ float lbessif_small(float const nu, float const z)
{
    USIZE_t M = 1000000, m;
    float x, f = 0.0, s = 0.0, so, y, g;
    x  = logf(0.5f*z);
    g  = lgammaf(nu+1.0f); /* gammaln(1) + gammaln(1+nu) */
    s  = 1.0f;
    f  = x * nu - g;
    for(m=1; m<=M; m++)
    {
        g += logf(m*(m + nu));
        y  = x * (2*m + nu) - g;
        if (y>f)
        {
            s *= expf(f - y);
            s += 1.0f;
            f  = y;
        }
        else
        {
            so = s;
            s += expf(y - f);
            if (so==s) break;
        }
        /* gammaln((m+1) + 1) + gammaln((m+1) + nu + 1) =
           gammaln( m    + 1) + gammaln( m    + nu + 1) + log(m + 1) + log(m + nu + 1) =
           gammaln( m    + 1) + gammaln( m    + nu + 1) + log((m + 1)*(m + nu + 1)) */
    }
    return (float)(logf(s) + f);
}


/* Compute log(bessi(nu, z))
   - requires z to be reasonably large for accurate results.
   - fails for nu0 = 0.0 (or close to 0.0).
*/
__device__ float lbessif_large(float const nu, float const z)
{
    static float tiny = 1.4014e-45f;
    float f, t, tt, ttt, us, den, tmp;
    f  = z/nu;
    f *= f;
    if (f>4.0f)
    {
        tmp = sqrtf(1.0f+1.0f/f);
        t   = z*tmp/nu;
        f   = nu*(t - logf(nu/z+tmp));
    }
    else
    {
        tmp = sqrtf(1.0f+f);
        t   = (tmp>1.0f) ? tmp : 1.0f;
        f   = nu*(t + logf(z/(nu*(1.0f+tmp))));
    }

    t    = 1.0f/t;
    tt   = t*t;
    ttt  = t*tt;
    us   = 1.0f;
    den  = nu;
    us  += t*(0.125f - tt*0.2083333333333333f)/den;
    den *= nu;
    us  += tt*(0.0703125f + tt*(-0.4010416666666667f + tt*0.3342013888888889f))/den;
    den *= nu;
    us  += ttt*(0.0732421875f + tt*(-0.8912109375f + tt*(1.846462673611111f - tt*1.025812596450617f)))/den;
    den *= nu;
    us  += tt*tt*(0.112152099609375f + tt*(-2.3640869140625f + tt*(8.78912353515625f +
           tt*(-11.20700261622299f + tt*4.669584423426248f))))/den;
    den *= nu;
    us  += tt*ttt*(0.2271080017089844f + tt*(-7.368794359479632f + tt*(42.53499874638846f +
           tt*(-91.81824154324002f + tt*(84.63621767460074f - tt*28.21207255820025f)))))/den;
    den *= nu;
    us  += ttt*ttt*(0.5725014209747314f + tt*(-26.49143048695155f + tt*(218.1905117442116f +
           tt*(-699.5796273761326f + tt*(1059.990452528f + tt*(-765.2524681411817f +
           tt*212.5701300392171f))))))/den;

    f   += 0.5*(logf((t < tiny) ? tiny : t) - logf(nu)) - 0.918938533204673f; /* 0.5*log(2*pi) */
    f   += logf((us<tiny) ? tiny : us);
    return f;
}


__device__ float lbessif(float const nu0, float const z0)
{
    float z = fabsf(z0), nu = fabsf(nu0);
    if (nu >= 15.0f)
        return lbessif_large(nu,z);
    else
    {
        float thr = 3.3333333333333333 * sqrtf(225.0f - nu*nu);
        if (z<thr)
            return lbessif_small(nu, z);
        else
        {
            nu = (nu < 1e-7) ? 1e-7 : nu;
            return lbessif_large(nu,z);
        }
    }
}
