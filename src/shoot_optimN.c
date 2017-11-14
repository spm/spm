/* $Id: shoot_optimN.c 7215 2017-11-14 15:57:55Z john $ */
/* (c) John Ashburner (2007) */

#include<mex.h>
#include<math.h>
extern double log(double x);
#define MAXD3 128

#include "shoot_boundary.h"
#include "shoot_multiscale.h"

static void choldc(int n, double a[], double p[])
{
    int i, j, k;
    double sm, sm0;

    sm0  = 1e-40;
    for(i=0; i<n; i++) sm0 = sm0 + a[i*n+i];
    sm0 *= 1e-7;
    sm0 *= sm0;
 /* for(i=0; i<n; i++) a[i*n+i] += sm0; */

    for(i=0; i<n; i++)
    {
        for(j=i; j<n; j++)
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

static void cholls(int n, double a[], double p[], double b[], double x[])
{
    int i, k;
    double sm;

    for(i=0; i<n; i++)
    {
        sm = b[i];
        for(k=i-1; k>=0; k--)
            sm -= a[i*n+k]*x[k];
        x[i] = sm/p[i];
    }
    for(i=n-1; i>=0; i--)
    {
        sm = x[i];
        for(k=i+1; k<n; k++)
            sm -= a[k*n+i]*x[k];
        x[i] = sm/p[i];
    }
}

static void Atimesp1(mwSize dm[], float A[], float p[], float Ap[])
{
    mwSize i,j, m = dm[0]*dm[1]*dm[2];
    float *pp[MAXD3], *pap[MAXD3], *pA[(MAXD3*(MAXD3+1))/2];

    for(i=0; i<dm[3]; i++)
    {
        pp[i]  = &p[m*i];
        pap[i] = &Ap[m*i];
    }
    for(i=0; i<(dm[3]*(dm[3]+1))/2; i++)
        pA[i] = &A[m*i];

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
                 pap[i][j] += ao*pp[k][j];
                 pap[k][j] += ao*pp[i][j];
             }
        }
    }
}

static void get_a(mwSize dm3, mwSize i, float *pa[],  double a[])
{
    mwSignedIndex m, n;
    mwSignedIndex o = dm3;
    for(m=0; m<dm3; m++)
    {
        a[m+dm3*m] = pa[m][i];
        for(n=m+1; n<dm3; n++,o++)
            a[n+dm3*m] = a[m+dm3*n] = pa[o][i];
    }
}

static double sumsq(mwSize dm[], float a[], float b[], double s[], double scal[], float u[])
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
    mwSignedIndex k;

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

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km2,km1,kp1,kp2;
        km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex i,m, jm2,jm1,jp1,jp2;
            float *p[MAXD3], *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
            double a1[MAXD3*MAXD3];

            for(m=0; m<dm[3]; m++)
            {
                pu[m]  = u+dm[0]*(j+dm[1]*(k+dm[2]*m));
                pb[m]  = b+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }

            if (a)
            {
                for(m=0; m<(dm[3]*(dm[3]+1))/2; m++)
                    pa[m]  = a+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }

            jm2 = (bound(j-2,dm[1])-j)*dm[0];
            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            jp2 = (bound(j+2,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex m, im2,im1,ip1,ip2;
                double tmp;

                im2 = bound(i-2,dm[0])-i;
                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;
                ip2 = bound(i+2,dm[0])-i;

                for(m=0; m<dm[3]; m++) p[m] = &(pu[m][i]);

                if (a) get_a(dm[3], i, pa, a1);

                for(m=0; m<dm[3]; m++)
                {
                    mwSignedIndex n;
                    float *pm =  p[m];
                    double pm0 = pm[0];
                    tmp =  (lam0*  pm0 +
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

                    if (a)
                    {
                        double *a11 = a1 + dm[3]*m;
                        for(n=0; n<dm[3]; n++) tmp += a11[n]*p[n][0];
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
    mwSignedIndex k;
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

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km2,km1,kp1,kp2;
        km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex i,m, jm2,jm1,jp1,jp2;
            float *pf[MAXD3], *pg[MAXD3];

            for(m=0; m<dm[3]; m++)
            {
                pf[m]  = f+dm[0]*(j+dm[1]*(k+dm[2]*m));
                pg[m]  = g+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }

            jm2 = (bound(j-2,dm[1])-j)*dm[0];
            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            jp2 = (bound(j+2,dm[1])-j)*dm[0];

            for(m=0; m<dm[3]; m++)
            {
                mwSignedIndex im2,im1,ip1,ip2;
                float *pf1 = pf[m], *pg1 = pg[m];
                for(i=0; i<dm[0]; i++)
                {
                    float *p = &pf1[i];
                    double p0 = p[0];

                    im2 = bound(i-2,dm[0])-i;
                    im1 = bound(i-1,dm[0])-i;
                    ip1 = bound(i+1,dm[0])-i;
                    ip2 = bound(i+2,dm[0])-i;
                    pg1[i] =(lam0*  p0 +
                           + w100*((p[im1        ]-p0) + (p[ip1        ]-p0))
                           + w010*((p[    jm1    ]-p0) + (p[    jp1    ]-p0))
                           + w001*((p[        km1]-p0) + (p[        kp1]-p0))
                           + w200*((p[im2        ]-p0) + (p[ip2        ]-p0))
                           + w020*((p[    jm2    ]-p0) + (p[    jp2    ]-p0))
                           + w002*((p[        km2]-p0) + (p[        kp2]-p0))
                           + w110*((p[im1+jm1    ]-p0) + (p[ip1+jm1    ]-p0) + (p[im1+jp1    ]-p0) + (p[ip1+jp1    ]-p0))
                           + w101*((p[im1    +km1]-p0) + (p[ip1    +km1]-p0) + (p[im1    +kp1]-p0) + (p[ip1    +kp1]-p0))
                           + w011*((p[    jm1+km1]-p0) + (p[    jp1+km1]-p0) + (p[    jm1+kp1]-p0) + (p[    jp1+kp1]-p0)))*scal[m];
                }
            }
        }
    }
}

static void relax(mwSize dm[], float a[], float b[], double s[], double scal[], int nit, float u[])
{
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

    w000 = w000*1.00001;

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
        mwSignedIndex i, j, k;
        for(k=(it/9)%3; k<dm[2]; k+=3)
        {
            mwSignedIndex km2, km1, kp1, kp2;
            km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
            kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];

            for(j=(it/3)%3; j<dm[1]; j+=3)
            {
                float *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
                double a1[MAXD3*MAXD3], cp[MAXD3], su[MAXD3];
                mwSignedIndex m, jm2,jm1,jp1,jp2;

                for(m=0; m<dm[3]; m++)
                {
                    pu[m]  = u+dm[0]*(j+dm[1]*(k+dm[2]*m));
                    pb[m]  = b+dm[0]*(j+dm[1]*(k+dm[2]*m));
                }

                if (a)
                {
                    for(m=0; m<(dm[3]*(dm[3]+1))/2; m++)
                        pa[m]  = a+dm[0]*(j+dm[1]*(k+dm[2]*m));
                }

                jm2 = (bound(j-2,dm[1])-j)*dm[0];
                jm1 = (bound(j-1,dm[1])-j)*dm[0];
                jp1 = (bound(j+1,dm[1])-j)*dm[0];
                jp2 = (bound(j+2,dm[1])-j)*dm[0];

                for(i=it%3; i<dm[0]; i+=3)
                {
                    mwSignedIndex im2,im1,ip1,ip2;

                    im2 = bound(i-2,dm[0])-i;
                    im1 = bound(i-1,dm[0])-i;
                    ip1 = bound(i+1,dm[0])-i;
                    ip2 = bound(i+2,dm[0])-i;

                    if (a) get_a(dm[3], i, pa, a1);

                    for(m=0; m<dm[3]; m++)
                    {
                        mwSignedIndex n;
                        float *pm  = &pu[m][i];
                        double pm0 = pm[0];
                        su[m] = pb[m][i]-
                               (lam0* pm0 
                              + w100*((pm[im1        ]-pm0) + (pm[ip1        ]-pm0))
                              + w010*((pm[    jm1    ]-pm0) + (pm[    jp1    ]-pm0))
                              + w001*((pm[        km1]-pm0) + (pm[        kp1]-pm0))
                              + w200*((pm[im2        ]-pm0) + (pm[ip2        ]-pm0))
                              + w020*((pm[    jm2    ]-pm0) + (pm[    jp2    ]-pm0))
                              + w002*((pm[        km2]-pm0) + (pm[        kp2]-pm0))
                              + w110*((pm[im1+jm1    ]-pm0) + (pm[ip1+jm1    ]-pm0) + (pm[im1+jp1    ]-pm0) + (pm[ip1+jp1    ]-pm0))
                              + w101*((pm[im1    +km1]-pm0) + (pm[ip1    +km1]-pm0) + (pm[im1    +kp1]-pm0) + (pm[ip1    +kp1]-pm0))
                              + w011*((pm[    jm1+km1]-pm0) + (pm[    jp1+km1]-pm0) + (pm[    jm1+kp1]-pm0) + (pm[    jp1+kp1]-pm0)))*scal[m];

                        if (a)
                            for(n=0; n<dm[3]; n++) su[m] -= a1[m*dm[3]+n]*pu[n][i];
                    }
                    if (a)
                    {
                        for(m=0; m<dm[3]; m++) a1[m+dm[3]*m] += w000*scal[m];
                        choldc(dm[3],a1,cp);
                        cholls(dm[3],a1,cp,su,su);
                        for(m=0; m<dm[3]; m++) pu[m][i] += su[m];
                    }
                    else
                    {
                        for(m=0; m<dm[3]; m++) pu[m][i] += su[m]/(w000*scal[m]);
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
    for(i=0; i<n; i++)
    {
        restrict_vol(na, a+i*na[0]*na[1]*na[2], nc, c+i*nc[0]*nc[1]*nc[2], b);
    }
}

static void prolong(mwSize n,  mwSize na[], float *a,  mwSize nc[], float *c, float *b)
{
    mwSignedIndex i;
    for(i=0; i<n; i++)
        resize_vol(na, a+i*na[0]*na[1]*na[2], nc, c+i*nc[0]*nc[1]*nc[2], b);
}

static void zeros(mwSize n, float *a)
{
    mwSignedIndex i;
    for(i=0; i<n; i++)
        a[i] = 0.0;
}

static void copy(mwSize n, float *a, float *b)
{
    mwSignedIndex i;
    for(i=0; i<n; i++)
        b[i] = a[i];
}

static void addto(mwSize n, float *a, float *b)
{
    mwSignedIndex i;
    for(i=0; i<n; i++)
        a[i] += b[i];
}

mwSignedIndex fmg_scratchsize(mwSize n0[])
{
    mwSignedIndex    n[32][3], m[32], bs, j;
    bs = 0;
    n[0][0] = n0[0];
    n[0][1] = n0[1];
    n[0][2] = n0[2];

    for(j=1; j<16; j++)
    {
        n[j][0] = ceil(n[j-1][0]/2.0);
        n[j][1] = ceil(n[j-1][1]/2.0);
        n[j][2] = ceil(n[j-1][2]/2.0);
        m[j]    = n[j][0]*n[j][1]*n[j][2];
        bs += m[j];
        if ((n[j][0]<2) && (n[j][1]<2) && (n[j][2]<2))
            break;
    }
    return((n0[3]*n0[0]*n0[1]*n0[2] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + (n0[3]*3+(n0[3]*(n0[3]+1))/2)*bs));
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
        n[j][0] = ceil(n[j-1][0]/2.0);
        n[j][1] = ceil(n[j-1][1]/2.0);
        n[j][2] = ceil(n[j-1][2]/2.0);
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
            int jj;
            for(jj=j; jj<ng-1; jj++)
            {
                relax(n[jj], a[jj], b[jj], param[jj], scal, nit, u[jj]);
                Atimesp(n[jj], a[jj], param[jj], scal, u[jj], res);
                for(i=0; i<n0[3]*m[jj]; i++)
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

