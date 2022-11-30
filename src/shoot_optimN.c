/*
 * John Ashburner
 * Copyright (C) 2011-2022 Wellcome Centre for Human Neuroimaging
 */

#include<math.h>
extern double log(double x);
#define MAXD3 128

#include "spm_mex.h"
#include "shoot_boundary.h"
#include "shoot_multiscale.h"

static void choldc(mwSize n, double a[], /*@out@*/ double p[])
{
    mwSignedIndex i, j, k;
    double sm, sm0;

    sm0  = 1e-40;
    for(i=0; i<(mwSignedIndex)n; i++) sm0 = sm0 + a[i*n+i];
    sm0 *= 1e-7;
    sm0 *= sm0;
 /* for(i=0; i<n; i++) a[i*n+i] += sm0; */

    for(i=0; i<(mwSignedIndex)n; i++)
    {
        for(j=i; j<(mwSignedIndex)n; j++)
        {
            sm = a[i*n+j];
            for(k=i-1; k>=0; k--)
               sm -= a[i*n+k] * a[j*n+k];
            if(i==j)
            {
                if(sm <= sm0) sm = sm0;
                p[i] = sqrt(sm);
            }
            else
                a[j*n+i] = sm / p[i];
        }
    }
}

static void cholls(mwSize n, double a[], double p[], /*@out@*/ double b[], /*@out@*/ double x[])
{
    mwSignedIndex i, k;
    double sm;

    for(i=0; i<(mwSignedIndex)n; i++)
    {
        sm = b[i];
        for(k=i-1; k>=0; k--)
            sm -= a[i*n+k]*x[k];
        x[i] = sm/p[i];
    }
    for(i=(mwSignedIndex)n-1; i>=0; i--)
    {
        sm = x[i];
        for(k=i+1; k<(mwSignedIndex)n; k++)
            sm -= a[k*n+i]*x[k];
        x[i] = sm/p[i];
    }
}

void Atimesp1(mwSize dm[], float A[], float p[], float Ap[])
{
    mwSize i, j, m = dm[0]*dm[1]*dm[2];
    float *pp[MAXD3], *pap[MAXD3], *pA[(MAXD3*(MAXD3+1))/2];

    for(i=0; i<dm[3]; i++)
    {
        pp[i]  = &p[m*i];
        pap[i] = &Ap[m*i];
    }
    for(i=0; i<(dm[3]*(dm[3]+1))/2; i++)
        pA[i] = &A[m*i];

#   pragma omp parallel for private(i)
    for(j=0; j<m; j++)
    {
        mwSize k, o;
        for(i=0; i<dm[3]; i++)
            pap[i][j] += pA[i][j]*pp[i][j];
        o = dm[3];
        for(i=0; i<dm[3]; i++)
        {
             for(k=i+1; k<dm[3]; k++,o++)
             {
                 double  ao = pA[o][j];
                 pap[i][j] += (float)(ao*pp[k][j]);
                 pap[k][j] += (float)(ao*pp[i][j]);
             }
        }
    }
}

static void get_a(mwSize dm3, mwSignedIndex i, /*@out@*/ float *pa[], /*@out@*/ double a[])
{
    mwSize m, n;
    mwSignedIndex o = (mwSignedIndex)dm3;
    for(m=0; m<dm3; m++)
    {
        a[m+dm3*m] = pa[m][i]*1.000001;
        for(n=m+1; n<dm3; n++,o++)
        {
            a[m+dm3*n] = pa[o][i];
            a[n+dm3*m] = pa[o][i];
        }
    }
}

/*@unused@*/ static double sumsq(mwSize dm[], float a[], float b[], double s[], double scal[], float u[])
{
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double lam0 = s[3], lam1 = s[4], lam2 = s[5];
    double ss = 0.0;
    mwSignedIndex j, k;

    w000 = lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2);
    w000 = w000 + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

#   pragma omp parallel for collapse(2) reduction(+: ss)
    for(k=0; k<(mwSignedIndex)dm[2]; k++)
    {
#       ifndef _OPENMP
            mwSignedIndex km2,km1,kp1,kp2;
            km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
            kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
#       endif

        for(j=0; j<(mwSignedIndex)dm[1]; j++)
        {
            mwSignedIndex i, m, jm2,jm1,jp1,jp2;
            float *p[MAXD3], *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
            double a1[MAXD3*MAXD3];
#           ifdef _OPENMP
                mwSignedIndex km2,km1,kp1,kp2;
                km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
                km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
                kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
                kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
#           endif

            for(m=0; m<(mwSignedIndex)dm[3]; m++)
            {
                pu[m]  = u+dm[0]*(j+dm[1]*(k+dm[2]*m));
                pb[m]  = b+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }

            if (a!=0)
            {
                for(m=0; m<(mwSignedIndex)(dm[3]*(dm[3]+1))/2; m++)
                    pa[m]  = a+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }

            jm2 = (bound(j-2,dm[1])-j)*dm[0];
            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            jp2 = (bound(j+2,dm[1])-j)*dm[0];

            for(i=0; i<(mwSignedIndex)dm[0]; i++)
            {
                mwSignedIndex im2,im1,ip1,ip2;
                double tmp;

                im2 = bound(i-2,dm[0])-i;
                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;
                ip2 = bound(i+2,dm[0])-i;

                for(m=0; m<(mwSignedIndex)dm[3]; m++) p[m] = &(pu[m][i]);

                if (a!=0) get_a(dm[3], i, pa, a1);

                for(m=0; m<(mwSignedIndex)dm[3]; m++)
                {
                    mwSignedIndex n;
                    float *pm =  p[m];
                    double pm0 = pm[0];
                    tmp =  (lam0*  pm0
                          + w100*((pm[im1        ]-pm0) + (pm[ip1        ]-pm0))
                          + w010*((pm[    jm1    ]-pm0) + (pm[    jp1    ]-pm0))
                          + w001*((pm[        km1]-pm0) + (pm[        kp1]-pm0))
                          + w200*((pm[im2        ]-pm0) + (pm[ip2        ]-pm0))
                          + w020*((pm[    jm2    ]-pm0) + (pm[    jp2    ]-pm0))
                          + w002*((pm[        km2]-pm0) + (pm[        kp2]-pm0))
                          + w110*((pm[im1+jm1    ]-pm0) + (pm[ip1+jm1    ]-pm0) + (pm[im1+jp1    ]-pm0) + (pm[ip1+jp1    ]-pm0))
                          + w101*((pm[im1    +km1]-pm0) + (pm[ip1    +km1]-pm0) + (pm[im1    +kp1]-pm0) + (pm[ip1    +kp1]-pm0))
                          + w011*((pm[    jm1+km1]-pm0) + (pm[    jp1+km1]-pm0) + (pm[    jm1+kp1]-pm0) + (pm[    jp1+kp1]-pm0)))*scal[m]
                          - pb[m][i];

/*
Note that there are numerical precision problems with this.
                    tmp =  (w000* pm[0] +
                          + w010*(pm[    jm1    ] + pm[    jp1    ])
                          + w020*(pm[    jm2    ] + pm[    jp2    ])
                          + w100*(pm[im1        ] + pm[ip1        ])
                          + w110*(pm[im1+jm1    ] + pm[ip1+jm1    ] + pm[im1+jp1    ] + pm[ip1+jp1    ])
                          + w200*(pm[im2        ] + pm[ip2        ])
                          + w001*(pm[        km1] + pm[        kp1])
                          + w101*(pm[im1    +km1] + pm[ip1    +km1] + pm[im1    +kp1] + pm[ip1    +kp1])
                          + w011*(pm[    jm1+km1] + pm[    jp1+km1] + pm[    jm1+kp1] + pm[    jp1+kp1])
                          + w002*(pm[        km2] + pm[        kp2]))*scal[m]
                          - pb[m][i];
*/

                    if (a!=0)
                    {
                        double *a11 = a1 + dm[3]*m;
                        for(n=0; n<(mwSignedIndex)dm[3]; n++) tmp += a11[n]*p[n][0];
                    }
                    ss += tmp*tmp;
                }
            }
        }
    }
    return(ss);
}

void LtLf(mwSize dm[], float f[], double s[], double scal[], float g[])
{
    mwSignedIndex j, k;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double lam0 = s[3], lam1 = s[4], lam2 = s[5];
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    w000 = lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2);
    w000 = w000 + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

    if (dm[0]<=2)
    {
        w000 += 2*w200;
        w200  = 0.0;
    }
    if (dm[1]<=2)
    {
        w000 += 2*w020;
        w020  = 0.0;
    }
    if (dm[2]<=2)
    {
        w000 += 2*w002;
        w002  = 0.0;
    }

    if (dm[0]==1)
    {
        w000 += 2*w100;
        w100  = 0.0;
        if (dm[1]==1)
        {
            w000 += 4*w110;
            w110  = 0.0;
        }
        if (dm[2]==1)
        {
            w000 += 4*w101;
            w101  = 0.0;
        }
    }
    if (dm[1]==1)
    {
        w000 += 2*w010;
        w010  = 0.0;
        if (dm[2]==1)
        {
            w000 += 4*w011;
            w011  = 0.0;
        }
    }
    if (dm[2]==1)
    {
        w000 += 2*w001;
        w001  = 0.0;
    }
    if (w000<0.0) w000=0.0;

#   pragma omp parallel for collapse(2)
    for(k=0; k<(mwSignedIndex)dm[2]; k++)
    {
#       ifndef _OPENMP
            mwSignedIndex km2,km1,kp1,kp2;
            km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
            kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
#       endif

        for(j=0; j<(mwSignedIndex)dm[1]; j++)
        {
            mwSignedIndex i,m,jm2,jm1,jp1,jp2;
            float *pf[MAXD3], *pg[MAXD3];
#           ifdef _OPENMP
                mwSignedIndex km2,km1,kp1,kp2;
                km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
                km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
                kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
                kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
#           endif

            for(m=0; m<(mwSignedIndex)dm[3]; m++)
            {
                pf[m]  = f+dm[0]*(j+dm[1]*(k+dm[2]*m));
                pg[m]  = g+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }

            jm2 = (bound(j-2,dm[1])-j)*dm[0];
            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            jp2 = (bound(j+2,dm[1])-j)*dm[0];

            for(m=0; m<(mwSignedIndex)dm[3]; m++)
            {
                float *pf1 = pf[m], *pg1 = pg[m];
                for(i=0; i<(mwSignedIndex)dm[0]; i++)
                {
                    float *p = &pf1[i];
                    double p0 = p[0];

                    mwSignedIndex im2,im1,ip1,ip2;

                    im2 = bound(i-2,dm[0])-i;
                    im1 = bound(i-1,dm[0])-i;
                    ip1 = bound(i+1,dm[0])-i;
                    ip2 = bound(i+2,dm[0])-i;
                    pg1[i] =(float)((lam0*  p0
                                   + w100*((p[im1        ]-p0) + (p[ip1        ]-p0))
                                   + w010*((p[    jm1    ]-p0) + (p[    jp1    ]-p0))
                                   + w001*((p[        km1]-p0) + (p[        kp1]-p0))
                                   + w200*((p[im2        ]-p0) + (p[ip2        ]-p0))
                                   + w020*((p[    jm2    ]-p0) + (p[    jp2    ]-p0))
                                   + w002*((p[        km2]-p0) + (p[        kp2]-p0))
                                   + w110*((p[im1+jm1    ]-p0) + (p[ip1+jm1    ]-p0) + (p[im1+jp1    ]-p0) + (p[ip1+jp1    ]-p0))
                                   + w101*((p[im1    +km1]-p0) + (p[ip1    +km1]-p0) + (p[im1    +kp1]-p0) + (p[ip1    +kp1]-p0))
                                   + w011*((p[    jm1+km1]-p0) + (p[    jp1+km1]-p0) + (p[    jm1+kp1]-p0) + (p[    jp1+kp1]-p0)))*scal[m]);
                }
            }
        }
    }
}

/* Non-stationary weighted version of LtLf.
 * It is used by TV/MTV regularisation, in a reweighted least-squares 
 * setting.
 * Currently, only membrane energy is implemented.
 *
 * . dm   {4}       - Volume dimensions (x/y/z/k)
 * . f    {dm(1:4)} - Input volume
 * . h    {dm(1:3)} - Weight volume
 * . s    {4}       - Regularisation parameters (vx/vy/vz/lam)
 * . scal {dm(4)}   - Regularisation scaling per feature
 * . g    {dm(1:4)} - Output voume
 */
void LtWLf(mwSize dm[], float f[], float h[], double s[], double scal[], float g[])
{
    mwSignedIndex i, j, k;
    double w000_000, w100_100, w010_010, w001_001;
    double lam = s[3];
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    /* Convolution kernels used to create final convolution weights.  */
    w000_000 =  lam*(v0 + v1 + v2);
    w100_100 = -lam*v0/2;
    w010_010 = -lam*v1/2;
    w001_001 = -lam*v2/2;

    /* Correct kernel if data is not 3D */
    if (dm[0]==1)
    {
        w000_000 -= lam*v0;
        w100_100  = 0.0;
        if (dm[1]==1)
        {
            w000_000 -= lam*v1;
            w010_010  = 0.0;
        }
        if (dm[2]==1)
        {
            w000_000 -= lam*v2;
            w001_001  = 0.0;
        }
    }
    else if (dm[1]==1)
    {
        w000_000 -= lam*v1;
        w010_010  = 0.0;
        if (dm[2]==1)
        {
            w000_000 -= lam*v2;
            w001_001  = 0.0;
        }
    }
    else if (dm[2]==1)
    {
        w000_000 -= lam*v2;
        w001_001  = 0.0;
    }
    if (w000_000<0.0) w000_000=0.0;
    
#   pragma omp parallel for collapse(3)
    for(k=0; k<(mwSignedIndex)dm[2]; k++)
    {
#       ifndef _OPENMP
            mwSignedIndex km1,kp1;
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
#       endif

        for(j=0; j<(mwSignedIndex)dm[1]; j++)
        {
#           ifndef _OPENMP
                mwSignedIndex jm1,jp1;
                jm1 = (bound(j-1,dm[1])-j)*dm[0];
                jp1 = (bound(j+1,dm[1])-j)*dm[0];
#           endif

            for(i=0; i<(mwSignedIndex)dm[0]; i++)
            {
                mwSignedIndex im1,ip1;
#               ifdef _OPENMP
                    mwSignedIndex jm1,jp1;
                    mwSignedIndex km1,kp1;
#               endif
                mwSignedIndex m;
                float *ph;
                double w1m00,w1p00,w01m0,w01p0,w001m,w001p;
                double ph0;

                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;

#               ifdef _OPENMP
                    jm1 = (bound(j-1,dm[1])-j)*dm[0];
                    jp1 = (bound(j+1,dm[1])-j)*dm[0];
                    km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
                    kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
#               endif

                /* Create convolution weights (that depend on h) */
                ph    = h+i+dm[0]*(j+dm[1]*k);
                ph0   = ph[0];
                w1m00 = w100_100*(ph0 + ph[im1]);
                w1p00 = w100_100*(ph0 + ph[ip1]);
                w01m0 = w010_010*(ph0 + ph[jm1]);
                w01p0 = w010_010*(ph0 + ph[jp1]);
                w001m = w001_001*(ph0 + ph[km1]);
                w001p = w001_001*(ph0 + ph[kp1]);


                for(m=0; m<(mwSignedIndex)dm[3]; m++)
                {
                    float *pf, *pg;
                    float pf0;
                    pf  = f+i+dm[0]*(j+dm[1]*(k+dm[2]*m));
                    pg  = g+i+dm[0]*(j+dm[1]*(k+dm[2]*m));
                    pf0 = pf[0];

                    /* Perform final convolution */
                    pg[0] = (float)(( w1m00*(pf[im1]-pf0) + w1p00*(pf[ip1]-pf0)
                                    + w01m0*(pf[jm1]-pf0) + w01p0*(pf[jp1]-pf0)
                                    + w001m*(pf[km1]-pf0) + w001p*(pf[kp1]-pf0))*scal[m]);
                }
            }
        }
    }
}

void solve(mwSize dm[], float a[], float b[], double s[], double scal[], float u[])
{
    double lam0 = s[3]; /* lam1 = s[4], lam2 = s[5]; */
    mwSignedIndex i, m;
    float *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
    double a1[MAXD3*MAXD3], cp[MAXD3], su[MAXD3];

    for(m=0; m<(mwSignedIndex)dm[3]; m++)
    {
        pu[m] = u+dm[0]*dm[1]*dm[2]*m;
        pb[m] = b+dm[0]*dm[1]*dm[2]*m;
    }
    if (a!=0)
    {
        for(m=0; m<(mwSignedIndex)(dm[3]*(dm[3]+1))/2; m++)
           pa[m] = a+dm[0]*dm[1]*dm[2]*m;
    }

#   pragma omp parallel for private(a1,cp,su)
    for(i=0; i<(mwSignedIndex)dm[0]*dm[1]*dm[2]; i++)
    {
        if (a!=0)
        {
            get_a(dm[3], i, pa, a1);
            for(m=0; m<(mwSignedIndex)dm[3]; m++)
            {
                su[m] = pb[m][i];
                a1[m+dm[3]*m] += lam0*scal[m];
            }
            choldc(dm[3],a1,cp);
            cholls(dm[3],a1,cp,su,su);
            for(m=0; m<(mwSignedIndex)dm[3]; m++) pu[m][i] = (float)(su[m]);
        }
        else
        {
            for(m=0; m<(mwSignedIndex)dm[3]; m++) pu[m][i] = (float)(pb[m][i]/(lam0*scal[m]));
        }
    }
}

/* Compute the diagonal of the inverse of the Hessian: diag(inv(H+L)).
 * To make the inversion tractable, L is approximated by its diagonal.
 * It is used for TV/MTV regularisation, in a (Bayesian) reweighted 
 * least-squares setting. Currently, only membrane energy is implemented.
 *
 * . dm   {4}         - Volume dimensions (x/y/z/f)
 * . a    {dm(1:3) k} - Hessian (k=f(f+1)/2)
 * . b    {dm(1:3)}   - Weight volume
 * . s    {4}         - Regularisation parameters (vx/vy/vz/lam)
 * . scal {f}         - Regularisation scaling per feature
 * . u    {dm(1:3) f} - Output volume
 *
 * If (u==0), the function returns the trace of the inverse of the 
 * Hessian: trace(inv(H+L)).
 */
double diaginv(mwSize dm[], float a[], float b[], double s[], double scal[], float u[])
{
    mwSignedIndex i, j, k;
    double w000_000, w000_100, w000_010, w000_001;
    double lam = s[3];
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double sum = 0;

    /* Convolution kernels used to create final convolution weights. */
    w000_000 = lam*(v0 + v1 + v2);
    w000_100 = lam*v0/2;
    w000_010 = lam*v1/2;
    w000_001 = lam*v2/2;

#   pragma omp parallel for collapse(3) reduction(+:sum)
    for(k=0; k<(mwSignedIndex)dm[2]; k++)
    {
#       ifndef _OPENMP
            mwSignedIndex km1,kp1;
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
#       endif

        for(j=0; j<(mwSignedIndex)dm[1]; j++)
        {
#           ifndef _OPENMP
                mwSignedIndex jm1,jp1;
                jm1 = (bound(j-1,dm[1])-j)*dm[0];
                jp1 = (bound(j+1,dm[1])-j)*dm[0];
#           endif

            for(i=0; i<(mwSignedIndex)dm[0]; i++)
            {
                mwSignedIndex im1,ip1;
#               ifdef _OPENMP
                    mwSignedIndex jm1,jp1;
                    mwSignedIndex km1,kp1;
#               endif
                mwSignedIndex m, mm;
                float *pb;
                double w000;
                double a1[MAXD3*MAXD3], cp[MAXD3], su[MAXD3];
                float *pa[(MAXD3*(MAXD3+1))/2];

                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;
#               ifdef _OPENMP
                    jm1 = (bound(j-1,dm[1])-j)*dm[0];
                    jp1 = (bound(j+1,dm[1])-j)*dm[0];
                    km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
                    kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
#               endif

                /* Create central convolution weight (that depends on h) */
                pb = b+i+dm[0]*(j+dm[1]*k);
                w000 =  w000_000*pb[0]
                     +  w000_100*(pb[im1] + pb[ip1])
                     +  w000_010*(pb[jm1] + pb[jp1])
                     +  w000_001*(pb[km1] + pb[kp1]);


                /* Add diagonal approximation of L to H */
                for(m=0; m<(mwSignedIndex)(dm[3]*(dm[3]+1))/2; m++)
                    pa[m] = a+dm[0]*(j+dm[1]*(k+dm[2]*m));
                get_a(dm[3], i, pa, a1);
                for(m=0; m<(mwSignedIndex)dm[3]; m++) a1[m+dm[3]*m] += w000*scal[m];
                /* Solve for inverse diagonal elements using Cholesky */
                choldc(dm[3],a1,cp);
                for(m=0; m<(mwSignedIndex)dm[3]; m++)
                {
                    for(mm=0; mm<(mwSignedIndex)dm[3]; mm++) su[mm] = 0;
                    su[m] = 1;
                    cholls(dm[3],a1,cp,su,su);
                    if(u!=0)
                    {
                        float *pu = u + i + dm[0]*(j+dm[1]*(k+dm[2]*m));
                        (*pu) = su[m];
                    }
                    else
                        sum += su[m];
                }

            }
        }
    }
    return(sum);
}

static void relax(mwSize dm[], float a[], float b[], double s[], double scal[], int nit, float u[])
{
    mwSignedIndex j, k;
    int it;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double lam0 = s[3], lam1 = s[4], lam2 = s[5];
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    w000 = (lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2));
    w000 = w000 + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

    w000 = w000*1.000001;

    if (dm[0]<=2)
    {
        w000 += 2*w200;
        w200  = 0.0;
    }
    if (dm[1]<=2)
    {
        w000 += 2*w020;
        w020  = 0.0;
    }
    if (dm[2]<=2)
    {
        w000 += 2*w002;
        w002  = 0.0;
    }

    if (dm[0]==1)
    {
        w000 += 2*w100;
        w100  = 0.0;
        if (dm[1]==1)
        {
            w000 += 4*w110;
            w110  = 0.0;
        }
        if (dm[2]==1)
        {
            w000 += 4*w101;
            w101  = 0.0;
        }
    }
    if (dm[1]==1)
    {
        w000 += 2*w010;
        w010  = 0.0;
        if (dm[2]==1)
        {
            w000 += 4*w011;
            w011  = 0.0;
        }
    }
    if (dm[2]==1)
    {
        w000 += 2*w001;
        w001  = 0.0;
    }
    if (w000<0.0) w000=0.0;

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("B%dx%dx%d: %g ", dm[0],dm[1],dm[2],sumsq(dm, a, b, s, scal, u));
#   endif

    for(it=0; it<27*nit; it++)
    {
#       pragma omp parallel for collapse(2)
        for(k=(it/9)%3; k<(mwSignedIndex)dm[2]; k+=3)
        {
#           ifndef _OPENMP
                mwSignedIndex km2, km1, kp1, kp2;
                km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
                km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
                kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
                kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
#           endif

            for(j=(it/3)%3; j<(mwSignedIndex)dm[1]; j+=3)
            {
                float *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
#               ifdef _OPENMP
                    mwSignedIndex km2, km1, kp1, kp2;
                    km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
                    km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
                    kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
                    kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
#               endif
                mwSignedIndex i, m, jm2,jm1,jp1,jp2;

                for(m=0; m<(mwSignedIndex)dm[3]; m++)
                {
                    pu[m]  = u+dm[0]*(j+dm[1]*(k+dm[2]*m));
                    pb[m]  = b+dm[0]*(j+dm[1]*(k+dm[2]*m));
                }

                if (a!=0)
                {
                    for(m=0; m<(mwSignedIndex)(dm[3]*(dm[3]+1))/2; m++)
                        pa[m]  = a+dm[0]*(j+dm[1]*(k+dm[2]*m));
                }

                jm2 = (bound(j-2,dm[1])-j)*dm[0];
                jm1 = (bound(j-1,dm[1])-j)*dm[0];
                jp1 = (bound(j+1,dm[1])-j)*dm[0];
                jp2 = (bound(j+2,dm[1])-j)*dm[0];

                for(i=it%3; i<(mwSignedIndex)dm[0]; i+=3)
                {
                    double a1[MAXD3*MAXD3], cp[MAXD3], su[MAXD3];
                    mwSignedIndex im2,im1,ip1,ip2;

                    im2 = bound(i-2,dm[0])-i;
                    im1 = bound(i-1,dm[0])-i;
                    ip1 = bound(i+1,dm[0])-i;
                    ip2 = bound(i+2,dm[0])-i;

                    if (a!=0) get_a(dm[3], i, pa, a1);
 
                    for(m=0; m<(mwSignedIndex)dm[3]; m++)
                    {
                        mwSignedIndex n;
                        float *pm  = &pu[m][i];
                        double pm0 = pm[0];
                        su[m] = (pb[m][i]-
                                       (lam0* pm0
                                      + w100*((pm[im1        ]-pm0) + (pm[ip1        ]-pm0))
                                      + w010*((pm[    jm1    ]-pm0) + (pm[    jp1    ]-pm0))
                                      + w001*((pm[        km1]-pm0) + (pm[        kp1]-pm0))
                                      + w200*((pm[im2        ]-pm0) + (pm[ip2        ]-pm0))
                                      + w020*((pm[    jm2    ]-pm0) + (pm[    jp2    ]-pm0))
                                      + w002*((pm[        km2]-pm0) + (pm[        kp2]-pm0))
                                      + w110*((pm[im1+jm1    ]-pm0) + (pm[ip1+jm1    ]-pm0) + (pm[im1+jp1    ]-pm0) + (pm[ip1+jp1    ]-pm0))
                                      + w101*((pm[im1    +km1]-pm0) + (pm[ip1    +km1]-pm0) + (pm[im1    +kp1]-pm0) + (pm[ip1    +kp1]-pm0))
                                      + w011*((pm[    jm1+km1]-pm0) + (pm[    jp1+km1]-pm0) + (pm[    jm1+kp1]-pm0) + (pm[    jp1+kp1]-pm0)))*scal[m]);

                        if (a!=0)
                            for(n=0; n<(mwSignedIndex)dm[3]; n++) su[m] -= a1[m*dm[3]+n]*pu[n][i];
                    }
                    if (a!=0)
                    {
                        for(m=0; m<(mwSignedIndex)dm[3]; m++) a1[m+dm[3]*m] += w000*scal[m];
                        choldc(dm[3],a1,cp);
                        cholls(dm[3],a1,cp,su,su);
                        for(m=0; m<(mwSignedIndex)dm[3]; m++)
                            pu[m][i] += su[m];
                    }
                    else
                    {
                        for(m=0; m<(mwSignedIndex)dm[3]; m++)
                            pu[m][i] += su[m]/(w000*scal[m]);
                    }
                }
            }
        }
#       ifdef VERBOSE
        if ((it%27) == 26)
            printf(" %g", sumsq(dm, a, b, s, scal, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}


static void Atimesp(mwSize dm[], float A[], double param[], double scal[], float p[], float Ap[])
{
    LtLf(dm, p, param, scal, Ap);
    Atimesp1(dm, A, p, Ap);
}


/*******************************************************/

static void restrictfcn(mwSize n,  mwSize na[], float *a,  mwSize nc[], float *c, float *b)
{
    mwSignedIndex i;
    for(i=0; i<(mwSignedIndex)n; i++)
    {
        restrict_vol(na, a+i*na[0]*na[1]*na[2], nc, c+i*nc[0]*nc[1]*nc[2], b);
    }
}

static void prolong(mwSize n,  mwSize na[], float *a,  mwSize nc[], float *c, float *b)
{
    mwSignedIndex i;
    for(i=0; i<(mwSignedIndex)n; i++)
        resize_vol(na, a+i*na[0]*na[1]*na[2], nc, c+i*nc[0]*nc[1]*nc[2], b);
}

static void zeros(mwSize n, float *a)
{
    mwSignedIndex i;
    for(i=0; i<(mwSignedIndex)n; i++)
        a[i] = 0.0;
}

static void copy(mwSize n, float *a, float *b)
{
    mwSignedIndex i;
    for(i=0; i<(mwSignedIndex)n; i++)
        b[i] = a[i];
}

static void addto(mwSize n, float *a, float *b)
{
    mwSignedIndex i;
    for(i=0; i<(mwSignedIndex)n; i++)
        a[i] += b[i];
}

mwSize fmg_scratchsize(mwSize n0[])
{
    mwIndex    n[32][3], m[32], bs, j;
    bs = 0;
    n[0][0] = n0[0];
    n[0][1] = n0[1];
    n[0][2] = n0[2];

    for(j=1; j<16; j++)
    {
        n[j][0] = (mwIndex)ceil((double)n[j-1][0]/2.0);
        n[j][1] = (mwIndex)ceil((double)n[j-1][1]/2.0);
        n[j][2] = (mwIndex)ceil((double)n[j-1][2]/2.0);
        m[j]    = n[j][0]*n[j][1]*n[j][2];
        bs += m[j];
        if ((n[j][0]<2) && (n[j][1]<2) && (n[j][2]<2))
            break;
    }
    return((mwSize)(n0[3]*n0[0]*n0[1]*n0[2] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + (n0[3]*3+(n0[3]*(n0[3]+1))/2)*bs));
}

/*
    Full Multigrid solver.  See Numerical Recipes (second edition) for more
    information
*/
void fmg(mwSize n0[], float *a0, float *b0, double param0[], double scal[], int c, int nit,
          float *u0, float *scratch)
{
    mwSignedIndex i, j, ng, bs;
    mwSize n[32][4], m[32];
    float *bo[32], *a[32], *b[32], *u[32], *res, *rbuf;
    double param[32][6];

#   ifdef VERBOSE
        printf("start=%g\n", sumsq(n0, a0, b0, param[0], scal, u0));
#   endif

    bo[0]   = b0;
    b[0]    = b0;
    u[0]    = u0;
    a[0]    = a0;
    n[0][0] = n0[0];
    n[0][1] = n0[1];
    n[0][2] = n0[2];
    n[0][3] = n0[3];
    m[0]    = n0[0]*n0[1]*n0[2];
    param[0][0] = param0[0];
    param[0][1] = param0[1];
    param[0][2] = param0[2];
    param[0][3] = param0[3];
    param[0][4] = param0[4];
    param[0][5] = param0[5];

    ng = 1;
    bs = 0;
    for(j=1; j<16; j++)
    {
        n[j][0] = (mwSize)ceil((double)n[j-1][0]/2.0);
        n[j][1] = (mwSize)ceil((double)n[j-1][1]/2.0);
        n[j][2] = (mwSize)ceil((double)n[j-1][2]/2.0);
        n[j][3] = n0[3];
        m[j]    = n[j][0]*n[j][1]*n[j][2];
        ng ++;
        bs += m[j];
        if ((n[j][0]<2) && (n[j][1]<2) && (n[j][2]<2))
            break;
    }

    res    = scratch;
    rbuf   = scratch + n0[3]*m[0];
    bo[1]  = scratch + n0[3]*m[0] + n[0][0]*n[1][1]+3*n0[0]*n0[1];
    b[1]   = scratch + n0[3]*m[0] + n[0][0]*n[1][1]+3*n0[0]*n0[1] + n0[3]*bs;
    u[1]   = scratch + n0[3]*m[0] + n[0][0]*n[1][1]+3*n0[0]*n0[1] + n0[3]*bs*2;
    a[1]   = scratch + n0[3]*m[0] + n[0][0]*n[1][1]+3*n0[0]*n0[1] + n0[3]*bs*3;

    for(j=2; j<ng; j++)
    {
        bo[j] = bo[j-1]+m[j-1]*n0[3];
        b[j]  =  b[j-1]+m[j-1]*n0[3];
        u[j]  =  u[j-1]+m[j-1]*n0[3];
        a[j]  =  a[j-1]+m[j-1]*(n0[3]*(n0[3]+1))/2;
    }
    for(j=1; j<ng; j++)
    {
        restrictfcn(n0[3],n[j-1],bo[j-1],n[j],bo[j],rbuf);
        restrictfcn((n0[3]*(n0[3]+1))/2,n[j-1],a[j-1],n[j],a[j],rbuf);

        param[j][0] = param0[0]*(double)n[j][0]/n0[0];
        param[j][1] = param0[1]*(double)n[j][1]/n0[1];
        param[j][2] = param0[2]*(double)n[j][2]/n0[2];
        param[j][3] = param[0][3];
        param[j][4] = param[0][4];
        param[j][5] = param[0][5];
    }
    relax(n[ng-1], a[ng-1], b[ng-1], param[ng-1], scal, nit, u[ng-1]);

    for(j=ng-2; j>=0; j--)
    {
        int jc;
        prolong(n0[3],n[j+1],u[j+1],n[j],u[j],rbuf);
        if(j>0) copy(n0[3]*m[j],bo[j],b[j]);
        for(jc=0; jc<c; jc++)
        {
            mwSignedIndex jj;
            for(jj=j; jj<ng-1; jj++)
            {
                relax(n[jj], a[jj], b[jj], param[jj], scal, nit, u[jj]);
                Atimesp(n[jj], a[jj], param[jj], scal, u[jj], res);
                for(i=0; i<(mwSignedIndex)n0[3]*m[jj]; i++)
                    res[i] = b[jj][i] - res[i];

                restrictfcn(n0[3],n[jj],res,n[jj+1],b[jj+1],rbuf);
                zeros(n0[3]*m[jj+1],u[jj+1]);
            }
            relax(n[ng-1], a[ng-1], b[ng-1], param[ng-1], scal, nit, u[ng-1]);

            for(jj=ng-2; jj>=j; jj--)
            {
                prolong(n0[3],n[jj+1],u[jj+1],n[jj],res,rbuf);
                addto(n0[3]*m[jj], u[jj], res);
                relax(n[jj], a[jj], b[jj], param[jj], scal, nit, u[jj]);
            }
        }
    }

/*  printf("end=%g\n", sumsq(n0, a0, b0, param[0], scal, u0)); */
}

/**************************************************************************
 *
 * SYMBOLIC DERIVATIONS FOR THE WEIGHTED CONVOLUTION PROBLEM (LtWLf)
 * (for MTV regularisation by reweighted least squares)
 *
 **************************************************************************

syms lam                % membrane penalty = 1/b^2 (b: Laplace parameter)
syms v0 v1 v2           % voxel size x/y/z (actually, 1/vx^2)
N  = 3;                 % Nb voxels in each dimension
% w  = sym(ones(N^3,1));  % This is the classical stationary case
w = sym('w', [N N N]);  % This is the new non-stationary case

O  = sym(zeros(N));
OO = kron(kron(O,O),O);
I  = sym(eye(N));
G1 = sym(spdiags(repmat([-1 1],N,1),[ 0 1],N,N)); G1(N,1) =  1; % forward difference
G2 = sym(spdiags(repmat([-1 1],N,1),[-1 0],N,N)); G2(1,N) = -1; % backward difference
G  = {G1 G2};

% Membrane energy
LL = sym(zeros(N^3,N^3));
for i=1:2
    Di = kron(I,kron(I,G{i}))*sqrt(v0); % 1st order / x
    LL = LL + lam*(Di.'*diag(w(:))*Di)/2;
end
for j=1:2
    Dj = kron(I,kron(G{j},I))*sqrt(v1); % 1st order / y
    LL = LL + lam*(Dj.'*diag(w(:))*Dj)/2;
end
for k=1:2
    Dk = kron(G{k},kron(I,I))*sqrt(v2); % 1st order / z
    LL = LL + lam*(Dk.'*diag(w(:))*Dk)/2;
end

% Reshape so that first 3: output voxels, second 3: input voxels
LL = reshape(LL, [N N N N N N]);
% Extract central output voxel
c = ceil(N/2);
LL = reshape(LL(c,c,c,:,:,:), [N N N]);
% The output value depends on all input neighbouring values
% We'll have some dependencies on the weight image 
LL = simplify(LL, 100);

% We have to construct 7 voxel-specific convolution weights.
% Each of these weights is obtained by convolving the weight image

wker = sym(zeros(N,N,N,N,N,N));
for i=1:N
for j=1:N
for k=1:N
    if ~isequal(LL(i,j,k), sym(0))
        for ii=1:N
        for jj=1:N
        for kk=1:N
            if ~isequal(LL(ii,jj,kk), sym(0))
                ww = sprintf('w%d_%d_%d', ii, jj, kk);
                wker(i,j,k,ii,jj,kk) = diff(LL(i,j,k),ww);
            end
        end
        end
        end
    end
end
end
end
wker = simplify(wker, 100);
c    = ceil(N/2);
k111 = reshape(wker(c,c,c,:,:,:), [N N N]);
k110 = reshape(wker(c,c,c-1,:,:,:), [N N N]);
k112 = reshape(wker(c,c,c+1,:,:,:), [N N N]);
k101 = reshape(wker(c,c-1,c,:,:,:), [N N N]);
k121 = reshape(wker(c,c+1,c,:,:,:), [N N N]);
k011 = reshape(wker(c-1,c,c,:,:,:), [N N N]);
k211 = reshape(wker(c+1,c,c,:,:,:), [N N N]);


% Convolution kernel for kernel weights:
%
% k000 (3x3x3)
% ------------
% kk000 = lam0 + 
%         lam1*(v0 + v1 + v2) + 
%         lam2*(4*v0^2 + 2*v0*v1 + 2*v0*v2 + 4*v1^2 + 2*v1*v2 + 4*v2^2)
% kk+00 = kk-00 = lam1*(v0/2) + lam2*(v0*v1 + v0*v2 + v0^2)
% kk0+0 = kk0-0 = lam1*(v1/2) + lam2*(v0*v1 + v1*v2 + v1^2)
% kk00+ = kk00- = lam1*(v2/2) + lam2*(v0*v2 + v1*v2 + v2^2)
% kk++0 = kk--0 = lam2*((v0*v1)/2)
% kk+0+ = kk-0- = lam2*((v0*v2)/2)
% kk0++ = kk0-- = lam2*((v1*v2)/2)
%
% 
% k+00 = sym(k-00) (3x3x3)
% ------------------------
% kk000 = kk+00 = lam1*(-v0/2) + lam2*(-(v0*(4*v0 + 2*v1 + 2*v2))/2)
% kk0+0 = kk0-0 = kk++0 = kk+-0 = lam2*(-(v0*v1)/2)
% kk00+ = kk00- = kk+0+ = kk+0- = lam2*(-(v0*v2)/2)
% kk-00 = 0
% kk-+0 = 0
% kk-0+ = 0
%
% k0+0 = sym(k0-0) (3x3x3)
% ------------------------
% kk000 = k0+0 = lam1*(-v1/2) + lam2*(-(v1*(2*v0 + 4*v1 + 2*v2))/2)
% kk+00 = kk-00 = kk++0 = kk-+0 = lam2*(-(v0*v1)/2)
% kk00+ = kk00- = kk0++ = kk0+- = lam2*(-(v1*v2)/2)
% kk0-0 = 0
% kk+-0 = 0
% kk0-+ = 0
%
% k00+ = sym(k00-) (3x3x3)
% ------------------------
% kk000 = k00+ = lam1*(-v2/2) + lam2*(-(v2*(2*v0 + 2*v1 + 4*v2))/2)
% kk00+ = kk00- = kk++0 = kk-+0 = lam2*(-(v0*v2)/2)
% kk00+ = kk00- = kk0++ = kk0+- = lam2*(-(v1*v2)/2)
% kk00- = 0
% kk+0- = 0
% kk0+- = 0
%
% k*00 = sym(k/00)
% ----------------
% kk+00 = lam2*v0^2
%
% k0*0 = sym(k0/0)
% ----------------
% kk0+0 = lam2*v1^2
%
% k00* = sym(k00/)
% ----------------
% kk00+ = lam2*v2^2
%
% k++0 = sym(k--0) = sym(k+-0) = sym(k-+0)
% ----------------------------------------
% kk000 = kk+00 = kk0+0 = kk++0 = lam2*((v0*v1)/2)
*/
