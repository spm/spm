/* $Id: shoot_dartel.c 7593 2019-05-20 18:58:16Z john $ */
/* (c) John Ashburner (2011) */

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "shoot_optim3d.h"
#include "shoot_diffeo3d.h"
#include "shoot_regularisers.h"
#include "shoot_boundary.h"

extern double   log(double x);
extern double   exp(double x);
#define LOG(x) (((x)>0) ? log(x+0.001): -6.9078)

/*
 * In place Cholesky decomposition
 */
void chol3(mwSize m, float A[])
{
    float *p00 = A,     *p11 = A+m,   *p22 = A+m*2,
          *p01 = A+m*3, *p02 = A+m*4, *p12 = A+m*5;
    double a00, a11, a22, a01, a02, a12;
    double s;
    mwSignedIndex i;
    for(i=0; i<m; i++)
    {
        a00 = *p00+1e-6;
        a11 = *p11+1e-6;
        a22 = *p22+1e-6;
        a01 = *p01;
        a02 = *p02;
        a12 = *p12;
        s        = sqrt(a00);
        *(p00++) = s;
        *(p01++) = a01/s;
        *(p02++) = a02/s;
        s        = a11 - a01*a01/a00;
        s        = sqrt(s);
        *(p11++) = s;
        s        = (a12 - a01*a02/a00)/s;
        *(p12++) = s;
        s        = a22 - a02*a02/a00 - s*s;
        *(p22++) = sqrt(s);
        /* printf("%g %g %g  %g %g %g\n", *(p00-1), *(p11-1), *(p22-1), *(p01-1), *(p02-1), *(p12-1)); */
    } 
}

/*
 * In place reconstruction from Cholesky decomposition
 */
void chol3recon(mwSize m, float A[])
{
    float *p00 = A,     *p11 = A+m,   *p22 = A+m*2,
          *p01 = A+m*3, *p02 = A+m*4, *p12 = A+m*5;
    double a00, a11, a22, a01, a02, a12;
    mwSignedIndex i;
    for(i=0; i<m; i++)
    {
        a00 = *p00;
        a11 = *p11;
        a22 = *p22;
        a01 = *p01;
        a02 = *p02;
        a12 = *p12;
        *(p00++) = a00*a00+a01*a01+a02*a02;
        *(p01++) = a00*a01+a01*a11+a02*a12;
        *(p02++) = a02*a00+a12*a01+a02*a22;
        *(p11++) = a01*a01+a11*a11+a12*a12;
        *(p12++) = a01*a02+a11*a12+a12*a22;
        *(p22++) = a02*a02+a12*a12+a22*a22;
    }
}

static mwSignedIndex pow2(int k)
{
    mwSignedIndex j0, td = 1;
    for(j0=0; j0<k; j0++)
        td = td*2;
    return(td);
}

/* Sample n points
 *  * s1 = f1(x,y,z)
 *   * s2 = f2(x,y,z)
 *    */
static void sampn_vox(mwSize dm[], float F[], mwSize n, mwSize mm, double x, double y, double z, double v[])
{
    mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
    mwSize j, o000, o100, o010, o110, o001, o101, o011, o111;
    mwSize tmpz, tmpy;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
    iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
    iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
    ix   = bound(ix  ,dm[0]);
    iy   = bound(iy  ,dm[1]);
    iz   = bound(iz  ,dm[2]);
    ix1  = bound(ix+1,dm[0]);
    iy1  = bound(iy+1,dm[1]);
    iz1  = bound(iz+1,dm[2]);

    tmpz  = dm[1]*iz;
    tmpy  = dm[0]*(iy + tmpz);
    o000  = (mwSize)(ix +tmpy);
    o100  = (mwSize)(ix1+tmpy);
    tmpy  = dm[0]*(iy1 + tmpz);
    o010  = (mwSize)(ix +tmpy);
    o110  = (mwSize)(ix1+tmpy);
    tmpz  = dm[1]*iz1;
    tmpy  = dm[0]*(iy + tmpz);
    o001  = (mwSize)(ix +tmpy);
    o101  = (mwSize)(ix1+tmpy);
    tmpy  = dm[0]*(iy1 + tmpz);
    o011  = (mwSize)(ix +tmpy);
    o111  = (mwSize)(ix1+tmpy);

    for(j=0; j<n; j++, F += mm)
    {
        v[j] = ((F[o000]*dx2 + F[o100]*dx1)*dy2 + (F[o010]*dx2 + F[o110]*dx1)*dy1)*dz2
             + ((F[o001]*dx2 + F[o101]*dx1)*dy2 + (F[o011]*dx2 + F[o111]*dx1)*dy1)*dz1;
    }
}

/*
 * t0 = Id + v0*sc
 * J0 = Id + |I+diag(v0)*sc|
 */
static void smalldef_jacdet(mwSize dm[], double sc, float v0[], float t0[], float J0[])
{
    mwSignedIndex j0, j1, j2;
    mwSignedIndex m = dm[0]*dm[1]*dm[2];
    double sc2 = sc/2.0;
    float *v1 = v0+m, *v2 = v1+m;
    
    for(j2=0; j2<dm[2]; j2++)
    {
        mwSignedIndex j2m1, j2p1;
        j2m1 = bound(j2-1,dm[2]);
        j2p1 = bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSignedIndex j1m1, j1p1;
            j1m1 = bound(j1-1,dm[1]);
            j1p1 = bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSignedIndex o, om1, op1;
                double j00,j10,j20, j01,j11,j21, j02,j12,j22;
                
                o         = j0+dm[0]*(j1+dm[1]*j2);
                t0[o    ] = (j0+1) + v0[o]*sc;
                t0[o+m  ] = (j1+1) + v1[o]*sc;
                t0[o+m*2] = (j2+1) + v2[o]*sc;

                om1 = bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                
                j00 = (v0[op1]-v0[om1])*sc2 + 1.0;
                j10 = (v1[op1]-v1[om1])*sc2;
                j20 = (v2[op1]-v2[om1])*sc2;

                om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                j01 = (v0[op1]-v0[om1])*sc2;
                j11 = (v1[op1]-v1[om1])*sc2 + 1.0;
                j21 = (v2[op1]-v2[om1])*sc2;

                om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                j02 = (v0[op1]-v0[om1])*sc2;
                j12 = (v1[op1]-v1[om1])*sc2;
                j22 = (v2[op1]-v2[om1])*sc2 + 1.0;

                J0[o] = j00*(j22*j11-j21*j12)
                      + j01*(j12*j20-j10*j22)
                      + j02*(j21*j10-j20*j11);
            }
        }
    }
}

/*
 * J0 := J0*inv(I+diag(v0)*sc)
 */
static void jac_div_smalldef(mwSize dm[], double sc, float v0[], float J0[])
{
    mwSignedIndex j0, j1, j2;
    mwSignedIndex m = dm[0]*dm[1]*dm[2];
    double sc2 = sc/2.0;
    float *v1 = v0+m, *v2 = v1+m;

    for(j2=0; j2<dm[2]; j2++)
    {
        mwSignedIndex j2m1, j2p1;
        j2m1 = bound(j2-1,dm[2]);
        j2p1 = bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSignedIndex j1m1, j1p1;
            j1m1 = bound(j1-1,dm[1]);
            j1p1 = bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSignedIndex o, om1, op1;
                double j00,j01,j02, j10,j11,j12, j20,j21,j22;
                double t00,t01,t02, t10,t11,t12, t20,t21,t22;
                double idt;

                om1 = bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                j00 = (v0[op1]-v0[om1])*sc2 + 1.0;
                j10 = (v1[op1]-v1[om1])*sc2;
                j20 = (v2[op1]-v2[om1])*sc2;

                om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                j01 = (v0[op1]-v0[om1])*sc2;
                j11 = (v1[op1]-v1[om1])*sc2 + 1.0;
                j21 = (v2[op1]-v2[om1])*sc2;

                om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                j02 = (v0[op1]-v0[om1])*sc2;
                j12 = (v1[op1]-v1[om1])*sc2;
                j22 = (v2[op1]-v2[om1])*sc2 + 1.0;

                /*
                syms j00 j01 j02 j10 j11 j12 j20 j21 j22
                syms d00 d01 d02 d10 d11 d12 d20 d21 d22
                J1 = [j00 j01 j02; j10 j11 j12; j20 j21 j22];
                inv(J1)
                J0 = [d00 d01 d02; d10 d11 d12; d20 d21 d22];
                J1*J0
                */

                t00 = j22*j11-j21*j12;
                t10 = j12*j20-j10*j22;
                t20 = j21*j10-j20*j11;
                t01 = j02*j21-j01*j22;
                t11 = j00*j22-j20*j02;
                t21 = j20*j01-j00*j21;
                t02 = j01*j12-j02*j11;
                t12 = j10*j02-j00*j12;
                t22 = j00*j11-j10*j01;
                idt = 1.0/(j00*t00+j01*t10+j02*t20);

                o   = j0+dm[0]*(j1+dm[1]*j2);
                j00 = J0[o    ]; j01 = J0[o+m*3]; j02 = J0[o+m*6];
                j10 = J0[o+m  ]; j11 = J0[o+m*4]; j12 = J0[o+m*7];
                j20 = J0[o+m*2]; j21 = J0[o+m*5]; j22 = J0[o+m*8];

                J0[o    ] = idt*(j00*t00+j01*t10+j02*t20);
                J0[o+m  ] = idt*(j10*t00+j11*t10+j12*t20);
                J0[o+m*2] = idt*(j20*t00+j21*t10+j22*t20);

                J0[o+m*3] = idt*(j00*t01+j01*t11+j02*t21);
                J0[o+m*4] = idt*(j10*t01+j11*t11+j12*t21);
                J0[o+m*5] = idt*(j20*t01+j21*t11+j22*t21);

                J0[o+m*6] = idt*(j00*t02+j01*t12+j02*t22);
                J0[o+m*7] = idt*(j10*t02+j11*t12+j12*t22);
                J0[o+m*8] = idt*(j20*t02+j21*t12+j22*t22);
            }
        }
    }
}

/*
 * Exponentiation with Jacobians
 */
void expdef(mwSize dm[], int k, double sc, float v[], float t0[], float t1[], float J0[], float J1[])
{
    float *optr;
    mwSignedIndex m = dm[0]*dm[1]*dm[2];
    int j;

    optr = t0;

    if(J0!=(float *)0)
    {
        smalldef_jac(dm, sc/pow2(k), v, t0, J0);
        for(j=0; j<k; j++)
        {
            float *tmpp;
            composition_jacobian(dm, dm[0]*dm[1]*dm[2], t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    }
    else
    {
        smalldef(dm, sc/pow2(k), v, t0);
        for(j=0; j<k; j++)
        {
            float *tmpp;
            composition(dm, dm[0]*dm[1]*dm[2], t0, t0, t1);
            tmpp = t0; t0   = t1; t1   = tmpp;
        }
    }
    if (optr != t0)
    {
        for(j=0; j<3*m; j++)
            t1[j] = t0[j];

        if (J0!=(float *)0)
            for(j=0; j<9*m; j++)
                J1[j] = J0[j];
    }
}

/*
 * Exponentiation with Jacobian determinants
 */
void expdefdet(mwSize dm[], int k, double sc, float v[], float t0[], float t1[], float J0[], float J1[])
{
    float *optr;
    mwSignedIndex m = dm[0]*dm[1]*dm[2];
    int j;

    optr = t0;

    if(J0!=(float *)0)
    {
        smalldef_jacdet(dm, sc/pow2(k), v, t0, J0);
        for(j=0; j<k; j++)
        {
            float *tmpp;
            composition_jacdet(dm, dm[0]*dm[1]*dm[2], t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    }
    else
    {
        smalldef(dm, sc/pow2(k), v, t0);
        for(j=0; j<k; j++)
        {
            float *tmpp;
            composition(dm, dm[0]*dm[1]*dm[2], t0, t0, t1);
            tmpp = t0; t0   = t1; t1   = tmpp;
        }
    }
    if (optr != t0)
    {
        for(j=0; j<3*m; j++)
            t1[j] = t0[j];

        if (J0!=(float *)0)
            for(j=0; j<m; j++)
                J1[j] = J0[j];
    }
}


static double smalldef_objfun_mn(mwSize dm[], float f[], float g[], float v[], float jd[], double sc, float b[], float A[])
{
    mwSignedIndex j, j0, j1, j2, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;

    j = 0;
    for(j2=0; j2<dm[2]; j2++)
        for(j1=0; j1<dm[1]; j1++)
            for(j0=0; j0<dm[0]; j0++, j++)
    {
        double x, y, z;
        mwSignedIndex ix, iy, iz, ix1, iy1, iz1, k;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx[128], dy[128], dz[128], Y[128], T[128], sT = 1.0, sY;
        double ta11, ta22, ta33, ta12, ta13, ta23;
        double tb1,  tb2,  tb3, tss;
        double sk000 = 1.0, sk100 = 1.0, sk010 = 1.0,
               sk110 = 1.0, sk001 = 1.0, sk101 = 1.0,
               sk011 = 1.0, sk111 = 1.0;
        
        x    = j0 + sc*v[j    ];
        y    = j1 + sc*v[j+m  ];
        z    = j2 + sc*v[j+m*2];
        ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = bound(ix,dm[0]);
        iy   = bound(iy,dm[1]);
        iz   = bound(iz,dm[2]);
        ix1  = bound(ix+1,dm[0]);
        iy1  = bound(iy+1,dm[1]);
        iz1  = bound(iz+1,dm[2]);
        sY   = 0.0;
        
        for(k=0; k<dm[3]; k++)
        {
            T[k]   = g[j + k*m];
            sT    -= T[k];
        }
        if (!mxIsFinite((double)sT))
        {
            A[j    ] = 0.0;
            A[j+m  ] = 0.0;
            A[j+m*2] = 0.0;
            A[j+m*3] = 0.0;
            A[j+m*4] = 0.0;
            A[j+m*5] = 0.0;
            b[j    ] = 0.0;
            b[j+m  ] = 0.0;
            b[j+m*2] = 0.0;
        }
        else
        {
            for(k=0; k<dm[3]; k++)
            {
                mwSignedIndex km = k*m;
                k000  = f[ix +dm[0]*(iy +dm[1]*iz ) + km];sk000-=k000;k000=LOG(k000);
                k100  = f[ix1+dm[0]*(iy +dm[1]*iz ) + km];sk100-=k100;k100=LOG(k100);
                k010  = f[ix +dm[0]*(iy1+dm[1]*iz ) + km];sk010-=k010;k010=LOG(k010);
                k110  = f[ix1+dm[0]*(iy1+dm[1]*iz ) + km];sk110-=k110;k110=LOG(k110);
                k001  = f[ix +dm[0]*(iy +dm[1]*iz1) + km];sk001-=k001;k001=LOG(k001);
                k101  = f[ix1+dm[0]*(iy +dm[1]*iz1) + km];sk101-=k101;k101=LOG(k101);
                k011  = f[ix +dm[0]*(iy1+dm[1]*iz1) + km];sk011-=k011;k011=LOG(k011);
                k111  = f[ix1+dm[0]*(iy1+dm[1]*iz1) + km];sk111-=k111;k111=LOG(k111);

                Y[k]  = exp(((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                          + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1);
                sY   += Y[k];

                dx[k] = -(((k000     - k100    )*dy2 + (k010     - k110    )*dy1)*dz2
                        + ((k001     - k101    )*dy2 + (k011     - k111    )*dy1)*dz1);
                dy[k] = -(((k000*dx2 + k100*dx1)     - (k010*dx2 + k110*dx1)    )*dz2
                        + ((k001*dx2 + k101*dx1)     - (k011*dx2 + k111*dx1)    )*dz1);
                dz[k] = -(((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)
                        - ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1));
            }

            k    = dm[3];
            T[k] = sT;
            k000 = LOG(sk000);
            k001 = LOG(sk001);
            k010 = LOG(sk010);
            k011 = LOG(sk011);
            k100 = LOG(sk100);
            k101 = LOG(sk101);
            k110 = LOG(sk110);
            k111 = LOG(sk111);

            Y[k] = exp(((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                     + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1);
            sY   += Y[k];

            dx[k] = -(((k000     - k100    )*dy2 + (k010     - k110    )*dy1)*dz2
                    + ((k001     - k101    )*dy2 + (k011     - k111    )*dy1)*dz1);
            dy[k] = -(((k000*dx2 + k100*dx1)     - (k010*dx2 + k110*dx1)    )*dz2
                    + ((k001*dx2 + k101*dx1)     - (k011*dx2 + k111*dx1)    )*dz1);
            dz[k] = -(((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)
                    - ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1));

            ta11 = ta22 = ta33 = ta12 = ta13 = ta23 = 0.0;
            tb1  = tb2  = tb3  = 0.0;
            tss  = 0.0;
            for(k=0; k<=dm[3]; k++)
            {
                Y[k] /= sY;
            }
            for(k=0; k<=dm[3]; k++)
            {
                double wt;
                mwSignedIndex k1;
                tss  += log(Y[k])*T[k];
                tb1  += (Y[k]-T[k])*dx[k];
                tb2  += (Y[k]-T[k])*dy[k];
                tb3  += (Y[k]-T[k])*dz[k];

                for(k1=0; k1<k; k1++)
                {
                    wt    =  -Y[k]*Y[k1];
                    ta11 += wt* dx[k]*dx[k1]*2;
                    ta22 += wt* dy[k]*dy[k1]*2;
                    ta33 += wt* dz[k]*dz[k1]*2;
                    ta12 += wt*(dx[k]*dy[k1] + dx[k1]*dy[k]);
                    ta13 += wt*(dx[k]*dz[k1] + dx[k1]*dz[k]);
                    ta23 += wt*(dy[k]*dz[k1] + dy[k1]*dz[k]);
                }
                wt    = Y[k]*(1.0-Y[k]);
                ta11 += wt*dx[k]*dx[k];
                ta22 += wt*dy[k]*dy[k];
                ta33 += wt*dz[k]*dz[k];
                ta12 += wt*dx[k]*dy[k];
                ta13 += wt*dx[k]*dz[k];
                ta23 += wt*dy[k]*dz[k];
            }

            if (jd != (float *)0)
            {
                double dt = jd[j];
                if (dt<0.0) dt = 0.0;
                A[j    ]  = ta11*dt;
                A[j+m  ]  = ta22*dt;
                A[j+m*2]  = ta33*dt;
                A[j+m*3]  = ta12*dt;
                A[j+m*4]  = ta13*dt;
                A[j+m*5]  = ta23*dt;
                b[j    ]  = tb1*dt;
                b[j+m  ]  = tb2*dt;
                b[j+m*2]  = tb3*dt;
                ssl      -= tss*dt;
            }
            else
            {
                A[j    ] = ta11;
                A[j+m  ] = ta22;
                A[j+m*2] = ta33;
                A[j+m*3] = ta12;
                A[j+m*4] = ta13;
                A[j+m*5] = ta23;
                b[j    ] = tb1;
                b[j+m  ] = tb2;
                b[j+m*2] = tb3;
                ssl     -= tss;
            }
        }
    }
    return(ssl);
}

static double smalldef_objfun2(mwSize dm[], float f[], float g[], float v[], float jd[], double sc, float b[], float A[])
{
    mwSignedIndex j, j0, j1, j2, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;

    j = 0;
    for(j2=0; j2<dm[2]; j2++)
        for(j1=0; j1<dm[1]; j1++)
            for(j0=0; j0<dm[0]; j0++, j++)
    {
        double x, y, z;
        mwSignedIndex ix, iy, iz, ix1, iy1, iz1, k;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        double d, dx, dy, dz, sd, sdx, sdy, sdz, ss;
        
        x    = j0 + sc*v[j    ];
        y    = j1 + sc*v[j+m  ];
        z    = j2 + sc*v[j+m*2];
        
        ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = bound(ix,dm[0]);
        iy   = bound(iy,dm[1]);
        iz   = bound(iz,dm[2]);
        ix1  = bound(ix+1,dm[0]);
        iy1  = bound(iy+1,dm[1]);
        iz1  = bound(iz+1,dm[2]);

        A[j    ] = 0.0;
        A[j+m  ] = 0.0;
        A[j+m*2] = 0.0;
        A[j+m*3] = 0.0;
        A[j+m*4] = 0.0;
        A[j+m*5] = 0.0;

        b[j    ] = 0.0;
        b[j+m  ] = 0.0;
        b[j+m*2] = 0.0;

        ss  = 0.0;
        sd  = 0.0;
        sdx = 0.0;
        sdy = 0.0;
        sdz = 0.0;

        for(k=0; k<dm[3]; k++)
        {
            mwSignedIndex km = k*m;
            k000  = f[ix +dm[0]*(iy +dm[1]*iz ) + km];
            k100  = f[ix1+dm[0]*(iy +dm[1]*iz ) + km];
            k010  = f[ix +dm[0]*(iy1+dm[1]*iz ) + km];
            k110  = f[ix1+dm[0]*(iy1+dm[1]*iz ) + km];
            k001  = f[ix +dm[0]*(iy +dm[1]*iz1) + km];
            k101  = f[ix1+dm[0]*(iy +dm[1]*iz1) + km];
            k011  = f[ix +dm[0]*(iy1+dm[1]*iz1) + km];
            k111  = f[ix1+dm[0]*(iy1+dm[1]*iz1) + km];

            d     = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                  + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 - g[j+km];

            dx    = ((k000     - k100    )*dy2 + (k010     - k110    )*dy1)*dz2
                  + ((k001     - k101    )*dy2 + (k011     - k111    )*dy1)*dz1;
            dy    = ((k000*dx2 + k100*dx1)     - (k010*dx2 + k110*dx1)    )*dz2
                  + ((k001*dx2 + k101*dx1)     - (k011*dx2 + k111*dx1)    )*dz1;
            dz    = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)
                  - ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1);

            sd  -= d;
            sdx -= dx;
            sdy -= dy;
            sdz -= dz;

            A[j    ] += dx*dx;
            A[j+m  ] += dy*dy;
            A[j+m*2] += dz*dz;
            A[j+m*3] += dx*dy;
            A[j+m*4] += dx*dz;
            A[j+m*5] += dy*dz;

            b[j    ] -= dx*d;
            b[j+m  ] -= dy*d;
            b[j+m*2] -= dz*d;

            ss += d*d;
        }
        A[j    ] += sdx*sdx;
        A[j+m  ] += sdy*sdy;
        A[j+m*2] += sdz*sdz;
        A[j+m*3] += sdx*sdy;
        A[j+m*4] += sdx*sdz;
        A[j+m*5] += sdy*sdz;

        b[j    ] -= sdx*sd;
        b[j+m  ] -= sdy*sd;
        b[j+m*2] -= sdz*sd;

        ss += sd*sd;
        
        if (jd != (float *)0)
        {
            double dt = jd[j];
            if (dt<0.0) dt = 0.0;
            A[j    ] *=dt;
            A[j+m  ] *=dt;
            A[j+m*2] *=dt;
            A[j+m*3] *=dt;
            A[j+m*4] *=dt;
            A[j+m*5] *=dt;
            b[j    ] *=dt;
            b[j+m  ] *=dt;
            b[j+m*2] *=dt;
            ss       *=dt;
        }
        ssl += ss;
    }
    return(0.5*ssl);
}

static double smalldef_objfun(mwSize dm[], float f[], float g[], float v[], float jd[], double sc, float b[], float A[])
{
    mwSignedIndex j,j0,j1,j2, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;

    if (dm[3]>1)
    {
        return(smalldef_objfun2(dm, f, g, v, jd, sc, b, A));
    }

    j = 0;
    for(j2=0; j2<dm[2]; j2++)
        for(j1=0; j1<dm[1]; j1++)
            for(j0=0; j0<dm[0]; j0++, j++)
    {
        double x, y, z;
        mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        double d, dx, dy, dz, dt = 1.0;
        
        x    = j0 + sc*v[j    ];
        y    = j1 + sc*v[j+m  ];
        z    = j2 + sc*v[j+m*2];

        ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = bound(ix,dm[0]);
        iy   = bound(iy,dm[1]);
        iz   = bound(iz,dm[2]);
        ix1  = bound(ix+1,dm[0]);
        iy1  = bound(iy+1,dm[1]);
        iz1  = bound(iz+1,dm[2]);
        
        k000  = f[ix +dm[0]*(iy +dm[1]*iz )];
        k100  = f[ix1+dm[0]*(iy +dm[1]*iz )];
        k010  = f[ix +dm[0]*(iy1+dm[1]*iz )];
        k110  = f[ix1+dm[0]*(iy1+dm[1]*iz )];
        k001  = f[ix +dm[0]*(iy +dm[1]*iz1)];
        k101  = f[ix1+dm[0]*(iy +dm[1]*iz1)];
        k011  = f[ix +dm[0]*(iy1+dm[1]*iz1)];
        k111  = f[ix1+dm[0]*(iy1+dm[1]*iz1)];
        
        d     = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
              + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 - g[j];
        dx    = ((k000     - k100    )*dy2 + (k010     - k110    )*dy1)*dz2
              + ((k001     - k101    )*dy2 + (k011     - k111    )*dy1)*dz1;
        dy    = ((k000*dx2 + k100*dx1)     - (k010*dx2 + k110*dx1)    )*dz2
              + ((k001*dx2 + k101*dx1)     - (k011*dx2 + k111*dx1)    )*dz1;
        dz    = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)
              - ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1);
        
        if (jd != (float *)0)
        {
            dt = jd[j];
            if (dt<0.0) dt = 0.0;
        }
        A[j    ] = dx*dx*dt;
        A[j+m  ] = dy*dy*dt;
        A[j+m*2] = dz*dz*dt;
        A[j+m*3] = dx*dy*dt;
        A[j+m*4] = dx*dz*dt;
        A[j+m*5] = dy*dz*dt;
        
        b[j    ] = -dx*d*dt;
        b[j+m  ] = -dy*d*dt;
        b[j+m*2] = -dz*d*dt;
        
        ssl += d*d*dt;
    }
    return(0.5*ssl);
}

static double initialise_objfun_mn(mwSize dm[], float f[], float g[], float t0[], float J0[], float jd[], float b[], float A[])
{
    mwSignedIndex j, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;
 
    for(j=0; j<m; j++)
    {
        double x, y, z;
        mwSignedIndex ix, iy, iz, ix1, iy1, iz1, k;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx0, dy0, dz0;
        double dx[128], dy[128], dz[128], Y[128], T[128], sT = 1.0, sY;
        double ta11, ta22, ta33, ta12, ta13, ta23;
        double tb1,  tb2,  tb3, tss;
        double sk000 = 1.0, sk100 = 1.0, sk010 = 1.0,
               sk110 = 1.0, sk001 = 1.0, sk101 = 1.0,
               sk011 = 1.0, sk111 = 1.0;

        x    = t0[j    ]-1.0;
        y    = t0[j+m  ]-1.0;
        z    = t0[j+m*2]-1.0;
        ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = bound(ix,dm[0]);
        iy   = bound(iy,dm[1]);
        iz   = bound(iz,dm[2]);
        ix1  = bound(ix+1,dm[0]);
        iy1  = bound(iy+1,dm[1]);
        iz1  = bound(iz+1,dm[2]);
        sY   = 0.0;
        
        for(k=0; k<dm[3]; k++)
        {
            T[k]   = g[j + k*m];
            sT    -= T[k];
        }
        if (!mxIsFinite((double)sT))
        {
            A[j    ] = 0.0;
            A[j+m  ] = 0.0;
            A[j+m*2] = 0.0;
            A[j+m*3] = 0.0;
            A[j+m*4] = 0.0;
            A[j+m*5] = 0.0;
            b[j    ] = 0.0;
            b[j+m  ] = 0.0;
            b[j+m*2] = 0.0;
        }
        else
        {
            for(k=0; k<dm[3]; k++)
            {
                mwSignedIndex km = k*m;
                k000  = f[ix +dm[0]*(iy +dm[1]*iz ) + km];sk000-=k000;k000=LOG(k000);
                k100  = f[ix1+dm[0]*(iy +dm[1]*iz ) + km];sk100-=k100;k100=LOG(k100);
                k010  = f[ix +dm[0]*(iy1+dm[1]*iz ) + km];sk010-=k010;k010=LOG(k010);
                k110  = f[ix1+dm[0]*(iy1+dm[1]*iz ) + km];sk110-=k110;k110=LOG(k110);
                k001  = f[ix +dm[0]*(iy +dm[1]*iz1) + km];sk001-=k001;k001=LOG(k001);
                k101  = f[ix1+dm[0]*(iy +dm[1]*iz1) + km];sk101-=k101;k101=LOG(k101);
                k011  = f[ix +dm[0]*(iy1+dm[1]*iz1) + km];sk011-=k011;k011=LOG(k011);
                k111  = f[ix1+dm[0]*(iy1+dm[1]*iz1) + km];sk111-=k111;k111=LOG(k111);

                Y[k]  = exp(((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                          + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1);
                sY   += Y[k];
            
                dx0   = ((k000     - k100    )*dy2 + (k010     - k110    )*dy1)*dz2
                      + ((k001     - k101    )*dy2 + (k011     - k111    )*dy1)*dz1;
                dy0   = ((k000*dx2 + k100*dx1)     - (k010*dx2 + k110*dx1)    )*dz2
                      + ((k001*dx2 + k101*dx1)     - (k011*dx2 + k111*dx1)    )*dz1;
                dz0   = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)
                      - ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1);

                dx[k] = -(J0[j    ]*dx0 + J0[j+  m]*dy0 + J0[j+2*m]*dz0);
                dy[k] = -(J0[j+3*m]*dx0 + J0[j+4*m]*dy0 + J0[j+5*m]*dz0);
                dz[k] = -(J0[j+6*m]*dx0 + J0[j+7*m]*dy0 + J0[j+8*m]*dz0);
            }

            k    = dm[3];
            T[k] = sT;

            k000 = LOG(sk000);
            k001 = LOG(sk001);
            k010 = LOG(sk010);
            k011 = LOG(sk011);
            k100 = LOG(sk100);
            k101 = LOG(sk101);
            k110 = LOG(sk110);
            k111 = LOG(sk111);

            Y[k]  = exp(((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                      + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1);
            sY   += Y[k];
            
            dx0   = ((k000     - k100    )*dy2 + (k010     - k110    )*dy1)*dz2
                  + ((k001     - k101    )*dy2 + (k011     - k111    )*dy1)*dz1;
            dy0   = ((k000*dx2 + k100*dx1)     - (k010*dx2 + k110*dx1)    )*dz2
                  + ((k001*dx2 + k101*dx1)     - (k011*dx2 + k111*dx1)    )*dz1;
            dz0   = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)
                  - ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1);

            dx[k] = -(J0[j    ]*dx0 + J0[j+  m]*dy0 + J0[j+2*m]*dz0);
            dy[k] = -(J0[j+3*m]*dx0 + J0[j+4*m]*dy0 + J0[j+5*m]*dz0);
            dz[k] = -(J0[j+6*m]*dx0 + J0[j+7*m]*dy0 + J0[j+8*m]*dz0);

            ta11 = ta22 = ta33 = ta12 = ta13 = ta23 = 0.0;
            tb1  = tb2  = tb3  = 0.0;
            tss  = 0.0;
            for(k=0; k<=dm[3]; k++)
            {
                double wt;
                mwSignedIndex k1;
                Y[k] /= sY;
                tss  += log(Y[k])*T[k];
                tb1  += (Y[k]-T[k])*dx[k];
                tb2  += (Y[k]-T[k])*dy[k];
                tb3  += (Y[k]-T[k])*dz[k];

                for(k1=0; k1<k; k1++)
                {
                    wt    =  -Y[k]*Y[k1];
                    ta11 += wt* dx[k]*dx[k1]*2;
                    ta22 += wt* dy[k]*dy[k1]*2;
                    ta33 += wt* dz[k]*dz[k1]*2;
                    ta12 += wt*(dx[k]*dy[k1] + dx[k1]*dy[k]);
                    ta13 += wt*(dx[k]*dz[k1] + dx[k1]*dz[k]);
                    ta23 += wt*(dy[k]*dz[k1] + dy[k1]*dz[k]);
                }
                wt    = Y[k]*(1.0-Y[k]);
                ta11 += wt*dx[k]*dx[k];
                ta22 += wt*dy[k]*dy[k];
                ta33 += wt*dz[k]*dz[k];
                ta12 += wt*dx[k]*dy[k];
                ta13 += wt*dx[k]*dz[k];
                ta23 += wt*dy[k]*dz[k];
            }
            if (jd != (float *)0)
            {
                double dt = jd[j];
                if (dt<0.0) dt = 0.0;
                A[j    ]  = ta11*dt;
                A[j+m  ]  = ta22*dt;
                A[j+m*2]  = ta33*dt;
                A[j+m*3]  = ta12*dt;
                A[j+m*4]  = ta13*dt;
                A[j+m*5]  = ta23*dt;
                b[j    ]  = tb1*dt;
                b[j+m  ]  = tb2*dt;
                b[j+m*2]  = tb3*dt;
                ssl      -= tss*dt;
            }
            else
            {
                A[j    ] = ta11;
                A[j+m  ] = ta22;
                A[j+m*2] = ta33;
                A[j+m*3] = ta12;
                A[j+m*4] = ta13;
                A[j+m*5] = ta23;
                b[j    ] = tb1;
                b[j+m  ] = tb2;
                b[j+m*2] = tb3;
                ssl     -= tss;
            }
        }
    }
    return(ssl);
}

static double initialise_objfun2(mwSize dm[], float f[], float g[], float t0[], float J0[], float jd[], float b[], float A[])
{
    mwSignedIndex j, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;
    
    for(j=0; j<m; j++)
    {
        double x, y, z;
        mwSignedIndex ix, iy, iz, ix1, iy1, iz1, k;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx0, dx1, dx2, dy0, dy1, dy2, dz0, dz1, dz2;
        double d, dx, dy, dz;
        double sd, sdx, sdy, sdz, ss;

        x    = t0[j    ]-1.0;
        y    = t0[j+m  ]-1.0;
        z    = t0[j+m*2]-1.0;
        ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = bound(ix,dm[0]);
        iy   = bound(iy,dm[1]);
        iz   = bound(iz,dm[2]);
        ix1  = bound(ix+1,dm[0]);
        iy1  = bound(iy+1,dm[1]);
        iz1  = bound(iz+1,dm[2]);

        A[j    ] = 0.0;
        A[j+m  ] = 0.0;
        A[j+m*2] = 0.0;
        A[j+m*3] = 0.0;
        A[j+m*4] = 0.0;
        A[j+m*5] = 0.0;

        b[j    ] = 0.0;
        b[j+m  ] = 0.0;
        b[j+m*2] = 0.0;

        ss  = 0.0;
        sd  = 0.0;
        sdx = 0.0;
        sdy = 0.0;
        sdz = 0.0;

        for(k=0; k<dm[3]; k++)
        {
            mwSignedIndex km = k*m;
            k000  = f[ix +dm[0]*(iy +dm[1]*iz ) + km];
            k100  = f[ix1+dm[0]*(iy +dm[1]*iz ) + km];
            k010  = f[ix +dm[0]*(iy1+dm[1]*iz ) + km];
            k110  = f[ix1+dm[0]*(iy1+dm[1]*iz ) + km];
            k001  = f[ix +dm[0]*(iy +dm[1]*iz1) + km];
            k101  = f[ix1+dm[0]*(iy +dm[1]*iz1) + km];
            k011  = f[ix +dm[0]*(iy1+dm[1]*iz1) + km];
            k111  = f[ix1+dm[0]*(iy1+dm[1]*iz1) + km];

            d     = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                  + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 - g[j+km];

            dx0   = ((k000     - k100    )*dy2 + (k010     - k110    )*dy1)*dz2
                  + ((k001     - k101    )*dy2 + (k011     - k111    )*dy1)*dz1;
            dy0   = ((k000*dx2 + k100*dx1)     - (k010*dx2 + k110*dx1)    )*dz2
                  + ((k001*dx2 + k101*dx1)     - (k011*dx2 + k111*dx1)    )*dz1;
            dz0   = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)
                  - ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1);

            dx   = -(J0[j    ]*dx0 + J0[j+  m]*dy0 + J0[j+2*m]*dz0);
            dy   = -(J0[j+3*m]*dx0 + J0[j+4*m]*dy0 + J0[j+5*m]*dz0);
            dz   = -(J0[j+6*m]*dx0 + J0[j+7*m]*dy0 + J0[j+8*m]*dz0);
            sd  -= d;
            sdx -= dx;
            sdy -= dy;
            sdz -= dz;

            A[j    ] += dx*dx;
            A[j+m  ] += dy*dy;
            A[j+m*2] += dz*dz;
            A[j+m*3] += dx*dy;
            A[j+m*4] += dx*dz;
            A[j+m*5] += dy*dz;

            b[j    ] += dx*d;
            b[j+m  ] += dy*d;
            b[j+m*2] += dz*d;

            ss       += d*d;
        }
        A[j    ] += sdx*sdx;
        A[j+m  ] += sdy*sdy;
        A[j+m*2] += sdz*sdz;
        A[j+m*3] += sdx*sdy;
        A[j+m*4] += sdx*sdz;
        A[j+m*5] += sdy*sdz;

        b[j    ] += sdx*sd;
        b[j+m  ] += sdy*sd;
        b[j+m*2] += sdz*sd;

        ss       += sd*sd;
        
        if (jd != (float *)0)
        {
            double dt = jd[j];
            if (dt<0.0) dt = 0.0;
            A[j    ] *=dt;
            A[j+m  ] *=dt;
            A[j+m*2] *=dt;
            A[j+m*3] *=dt;
            A[j+m*4] *=dt;
            A[j+m*5] *=dt;
            b[j    ] *=dt;
            b[j+m  ] *=dt;
            b[j+m*2] *=dt;
            ss       *=dt;
        }
        ssl += ss;
    }
    return(0.5*ssl);
}

static double initialise_objfun(mwSize dm[], float f[], float g[], float t0[], float J0[], float jd[], float b[], float A[])
{
    mwSignedIndex j, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0, dt = 1.0;

    if (dm[3]>1)
    {
        return(initialise_objfun2(dm, f, g, t0, J0, jd, b, A));
    }

    for(j=0; j<m; j++)
    {
        double x, y, z;
        mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx0, dx1, dx2, dy0, dy1, dy2, dz0, dz1, dz2;
        double d, dx, dy, dz;

        x    = t0[j    ]-1.0;
        y    = t0[j+m  ]-1.0;
        z    = t0[j+m*2]-1.0;
        ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = bound(ix,dm[0]);
        iy   = bound(iy,dm[1]);
        iz   = bound(iz,dm[2]);
        ix1  = bound(ix+1,dm[0]);
        iy1  = bound(iy+1,dm[1]);
        iz1  = bound(iz+1,dm[2]);

        k000  = f[ix +dm[0]*(iy +dm[1]*iz )];
        k100  = f[ix1+dm[0]*(iy +dm[1]*iz )];
        k010  = f[ix +dm[0]*(iy1+dm[1]*iz )];
        k110  = f[ix1+dm[0]*(iy1+dm[1]*iz )];
        k001  = f[ix +dm[0]*(iy +dm[1]*iz1)];
        k101  = f[ix1+dm[0]*(iy +dm[1]*iz1)];
        k011  = f[ix +dm[0]*(iy1+dm[1]*iz1)];
        k111  = f[ix1+dm[0]*(iy1+dm[1]*iz1)];

        d     = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
              + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 - g[j];
        dx0   = ((k000     - k100    )*dy2 + (k010     - k110    )*dy1)*dz2
              + ((k001     - k101    )*dy2 + (k011     - k111    )*dy1)*dz1;
        dy0   = ((k000*dx2 + k100*dx1)     - (k010*dx2 + k110*dx1)    )*dz2
              + ((k001*dx2 + k101*dx1)     - (k011*dx2 + k111*dx1)    )*dz1;
        dz0   = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)
              - ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1);

        dx   = -(J0[j    ]*dx0 + J0[j+  m]*dy0 + J0[j+2*m]*dz0);
        dy   = -(J0[j+3*m]*dx0 + J0[j+4*m]*dy0 + J0[j+5*m]*dz0);
        dz   = -(J0[j+6*m]*dx0 + J0[j+7*m]*dy0 + J0[j+8*m]*dz0);

        if (jd != (float *)0)
        {
            dt = jd[j];
            if (dt<0.0) dt = 0.0;
        }

        A[j    ] = dx*dx*dt;
        A[j+m  ] = dy*dy*dt;
        A[j+m*2] = dz*dz*dt;
        A[j+m*3] = dx*dy*dt;
        A[j+m*4] = dx*dz*dt;
        A[j+m*5] = dy*dz*dt;

        b[j    ] = dx*d*dt;
        b[j+m  ] = dy*d*dt;
        b[j+m*2] = dz*d*dt;

        ssl += d*d*dt;
    }
    return(0.5*ssl);
}


static void squaring(mwSize dm[], int k, int save_transf, float b[], float A[], float t0[], float t1[], float J0[], float J1[])
{
    mwSignedIndex i, j, m = dm[0]*dm[1]*dm[2];
    float *ptr = t0;

    for(i=0; i<k; i++)
    {
        float *buf1, *buf2;
        buf1 = t1; /* Re-use some memory */
        buf2 = J1;

#ifdef CHOL
        chol3(m, A);
#endif

        for(j=0; j<m; j++)
        {
            double x, y, z;
            double j00, j01, j02, j10, j11, j12, j20, j21, j22, dt;
            double a00, a11, a22, a01, a02, a12;
            double b0, b1, b2, tmp0, tmp1, tmp2;
            double as[6], bs[3];

            /*
            syms j00 j01 j02 j10 j11 j12 j20 j21 j22
            syms a00 a11 a22 a01 a02 a12
            syms b0 b1 b2
            J = [j00 j01 j02; j10 j11 j12; j20 j21 j22];
            A = [a00 a01 a02; a01 a11 a12; a02 a12 a22];
            b = [b0; b1; b2];
            J.'*b
            J.'*A*J
            */

            x   = t0[j    ]-1.0;
            y   = t0[j+m  ]-1.0;
            z   = t0[j+m*2]-1.0;

            j00 = J0[j    ]; j01 = J0[j+m*3]; j02 = J0[j+m*6];
            j10 = J0[j+m  ]; j11 = J0[j+m*4]; j12 = J0[j+m*7];
            j20 = J0[j+m*2]; j21 = J0[j+m*5]; j22 = J0[j+m*8];

            dt  = j00*(j11*j22-j12*j21)+j10*(j02*j21-j01*j22)+j20*(j01*j12-j02*j11);

            /* J'*b */
            sampn_vox(dm, b, 3, m, x, y, z, bs);
            b0 = bs[0];
            b1 = bs[1];
            b2 = bs[2];

            buf1[j    ] = dt*(b0*j00+b1*j10+b2*j20);
            buf1[j+m  ] = dt*(b0*j01+b1*j11+b2*j21);
            buf1[j+m*2] = dt*(b0*j02+b1*j12+b2*j22);

            /* J'*A*J */
            sampn_vox(dm, A, 6, m, x, y, z, as);
            a00 = as[0];
            a11 = as[1];
            a22 = as[2];
            a01 = as[3];
            a02 = as[4];
            a12 = as[5];

            /* rearranged for speed */
            tmp0        = j00*a00+j10*a01+j20*a02;
            tmp1        = j00*a01+j10*a11+j20*a12;
            tmp2        = j00*a02+j10*a12+j20*a22;
            buf2[j    ] = dt*(tmp0*j00+tmp1*j10+tmp2*j20);
            buf2[j+m*3] = dt*(tmp0*j01+tmp1*j11+tmp2*j21);
            buf2[j+m*4] = dt*(tmp0*j02+tmp1*j12+tmp2*j22);

            tmp0        = j01*a00+j11*a01+j21*a02;
            tmp1        = j01*a01+j11*a11+j21*a12;
            tmp2        = j01*a02+j11*a12+j21*a22;
            buf2[j+m  ] = dt*(tmp0*j01+tmp1*j11+tmp2*j21);
            buf2[j+m*5] = dt*(tmp0*j02+tmp1*j12+tmp2*j22);

            buf2[j+m*2] = dt*((j02*a00+j12*a01+j22*a02)*j02+(j02*a01+j12*a11+j22*a12)*j12+(j02*a02+j12*a12+j22*a22)*j22);
        }

#ifdef CHOL
        chol3recon(m, A);
        chol3recon(m, buf2);
#endif

        for(j=0; j<m*3; j++) b[j] += buf1[j];
        for(j=0; j<m*6; j++) A[j] += buf2[j];
        if (save_transf || (i<k-1))
        {
            float *tmpp;
            composition_jacobian(dm, dm[0]*dm[1]*dm[2], t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    }
    if (save_transf && ptr!=t0)
    {
        for(j=0; j<m*3; j++)
            t1[j] = t0[j];
        for(j=0; j<m*9; j++)
            J1[j] = J0[j];
    }
}

mwSize iteration_scratchsize(mwSize dm[], int code, int k)
{
    mwSignedIndex m1, m2;
    mwSignedIndex m = dm[0]*dm[1]*dm[2];
    if (k>0)
    {
        m1 = 30*m;
        if (code==1) m1 += 9*m;
        m2 = 9*m+fmg3_scratchsize(dm,1);
        if (m1>m2)
            return(m1);
        else
            return(m2);
    }
    else
    {
        m1 = 9*m;
        if (code==1) m1 += 6*m;
        m2 = 9*m + fmg3_scratchsize(dm,1);
        if (m1>m2)
            return(m1);
        else
            return(m2);
    }
}

void iteration(mwSize dm[], int k, float v[], float g[], float f[], float jd[],
               double param0[], double lmreg0, int cycles, int its, int code,
               float ov[], double ll[], float *buf)
{
    float *sbuf;
    float *b, *A, sc;
    double ssl, ssp;
    static double param[] = {1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0};
    mwSignedIndex m = dm[0]*dm[1]*dm[2];
    mwSignedIndex j;

    /*
        Allocate memory.
          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
        [ A  A  A  A  A  A  t  t  t  J  J  J  J  J  J  J  J  J  t  t  t  J  J  J  J  J  J  J  J  J] for computing derivatives
    */
    b    = ov;
    A    = buf;
    sbuf = buf +  6*m;

    if(k>0)
    {
        float *t0, *t1, *J0, *J1;
        t0   = buf +  6*m;
        J0   = buf +  9*m;
        t1   = buf + 18*m;
        J1   = buf + 21*m;

        sc = 1.0/pow2(k);
        expdef(dm, k, 1.0, v, t0, t1, J0, J1);
        jac_div_smalldef(dm, sc, v, J0);
        if (code==2)
            ssl = initialise_objfun_mn(dm, f, g, t0, J0, jd, b, A);
        else
            ssl = initialise_objfun(dm, f, g, t0, J0, jd, b, A);
        smalldef_jac(dm, -sc, v, t0, J0);
        squaring(dm, k, code==1, b, A, t0, t1, J0, J1);
        if (code==1)
        {
            float *b1, *A1;
            A1   = buf + 30*m;
            b1   = buf + 36*m;
            jac_div_smalldef(dm, -sc, v, J0);
            ssl += initialise_objfun(dm, g, f, t0, J0, (float *)0, b1, A1);
            smalldef_jac(dm, sc, v, t0, J0);
            squaring(dm, k, 0, b1, A1, t0, t1, J0, J1);
            for(j=0; j<m*3; j++) b[j] -= b1[j];
            for(j=0; j<m*6; j++) A[j] += A1[j];
        }
    }
    else
    {
        sc  = 1.0;
        if (code==2)
            ssl = smalldef_objfun_mn(dm, f, g, v, jd, 1.0, b, A);
        else
            ssl = smalldef_objfun(dm, f, g, v, jd, 1.0, b, A);
        if (code==1)
        {
            float *b1, *A1;
            A1   = buf + 6*m;
            b1   = buf + 12*m;
            ssl += smalldef_objfun(dm, g, f, v, (float *)0, -1.0, b1, A1);
            for(j=0; j<m*3; j++) b[j] -= b1[j];
            for(j=0; j<m*6; j++) A[j] += A1[j];
        }
    }

    param[3] = param0[3];
    param[4] = param0[4];
    param[5] = param0[5];
    param[6] = param0[6];
    param[7] = param0[7];

    vel2mom(dm, v, param, sbuf);

    ssp = 0.0;
    for(j=0; j<m*3; j++)
    {
        b[j] = b[j]*sc + sbuf[j];
        ssp += sbuf[j]*v[j];
    }

    ll[0] = ssl;
    ll[1] = ssp*0.5;
    ll[2] = norm(m*3,b);

    for(j=0; j<m*6; j++) A[j] *= sc;

    /* Solve equations for Levenberg-Marquardt update:
     * v = v - inv(H + L'*L + R)*(d + L'*L*v)
     *     v: velocity or flow field
     *     H: matrix of second derivatives
     *     L: regularisation (L'*L is the inverse of the prior covariance)
     *     R: Levenberg-Marquardt regularisation
     *     d: vector of first derivatives
     */

    if (lmreg0>0.0) param[3] = param[3] + lmreg0;

    fmg3(dm, A, b, param, cycles, its, sbuf, sbuf+3*m); 
    for(j=0; j<m*3; j++) ov[j] = v[j] - sbuf[j];
}

/*
 * Attempt to unwrap the deformations.
 * Note: this is not always guaranteed to work,
 * but it should for most cases.
 */
static void unwrap(mwSize dm[], float f[])
{
    mwSize i0, i1, i2;

    if (get_bound()!=0)
        return;

    for(i2=0; i2<dm[2]; i2++)
    {
        float *pt = f + (i2+2*dm[2])*dm[0]*dm[1];
        if (i2==0)
        {
            for(i1=0; i1<dm[1]*dm[0]; i1++)
                pt[i1] = pt[i1]-(float)floor((double)(pt[i1]/dm[2])+0.5)*(float)dm[2];
        }
        else
        {
            for(i1=0; i1<dm[1]*dm[0]; i1++)
                pt[i1] = pt[i1]-(float)floor((double)((pt[i1]-pt[i1-dm[0]*dm[1]])/dm[2])+0.5)*(float)dm[2];
        }
    }

    for(i1=0; i1<dm[1]; i1++)
    {
        float *pt = f + (i1+dm[2]*dm[1])*dm[0];
        if (i1==0)
        {
            for(i2=0; i2<dm[2]; i2++)
            {
                float *pPsi1 = pt+i2*dm[0]*dm[1];
                for(i0=0; i0<dm[0]; i0++)
                {
                    pPsi1[i0] = pPsi1[i0]-(float)floor((double)(pPsi1[i0]/dm[1])+0.5)*(float)dm[1];
                }
            }
        }
        else
        {
            for(i2=0; i2<dm[2]; i2++)
            {
                float *pPsi1 = pt+i2*dm[0]*dm[1];
                for(i0=0; i0<dm[0]; i0++)
                {
                    pPsi1[i0] = pPsi1[i0]-(float)floor((double)((pPsi1[i0]-pPsi1[i0-dm[0]])/dm[1])+0.5)*(float)dm[1];
                }
            }
        }
    }

    for(i0=0; i0<dm[0]; i0++)
    {
        float *pt = f+i0;
        if (i0==0)
        {
            for(i2=0; i2<dm[2]; i2++)
            {
                float *pPsi1 = pt + i2*dm[0]*dm[1];
                for(i1=0; i1<dm[0]*dm[1]; i1+=dm[0])
                    pPsi1[i1] = pPsi1[i1]-(float)floor((double)(pPsi1[i1]/dm[0])+0.5)*(float)dm[0];
            }
        }
        else
        {
            for(i2=0; i2<dm[2]; i2++)
            {
                float *pPsi1 = pt + i2*dm[0]*dm[1];
                for(i1=0; i1<dm[0]*dm[1]; i1+=dm[0])
                    pPsi1[i1] = pPsi1[i1]-(float)floor((double)((pPsi1[i1]-pPsi1[i1-1])/dm[0])+0.5)*(float)dm[0];
            }
        }
    }
}


void dartel_mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[])
{
    int        i, k=10, cycles=4, its=2, code=0;
    mwSize     dm[5];
    double     lmreg0=0.0, *ll;
    float      *v, *g, *f, *jd = (float *)0, *ov, *scratch;
    static double param[] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    static mwSize nll[] = {1, 3, 1};

    if ((nrhs!=4 && nrhs!=5) || nlhs>2)
        mexErrMsgTxt("Incorrect usage");

    for(i=0; i<3; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
            mexErrMsgTxt("Data must be numeric, real, full and double");

    if (mxGetNumberOfDimensions(prhs[0])!=4) mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetNumberOfDimensions(prhs[1])>4) mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetNumberOfDimensions(prhs[2])!=mxGetNumberOfDimensions(prhs[1])) mexErrMsgTxt("Incompatible number of dimensions.");
    dm[0] = mxGetDimensions(prhs[0])[0];
    dm[1] = mxGetDimensions(prhs[0])[1];
    dm[2] = mxGetDimensions(prhs[0])[2];
    dm[3] = mxGetDimensions(prhs[0])[3];

    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension of 1st arg must be 3.");

    if (mxGetDimensions(prhs[1])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[1])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");
    if (mxGetNumberOfDimensions(prhs[1])>=3 && mxGetDimensions(prhs[1])[2] != dm[2])
        mexErrMsgTxt("Incompatible 3rd dimension.");

    if (mxGetDimensions(prhs[2])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[2])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");
    if (mxGetNumberOfDimensions(prhs[2])>=3 && mxGetDimensions(prhs[2])[2] != dm[2])
        mexErrMsgTxt("Incompatible 3rd dimension.");

    if (nrhs>=5)
    {
        if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) || mxIsSparse(prhs[4]) || !mxIsSingle(prhs[4]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
        if (mxGetNumberOfDimensions(prhs[4])!=3) mexErrMsgTxt("Wrong number of dimensions.");
        if (mxGetDimensions(prhs[4])[0] != dm[0])
            mexErrMsgTxt("Incompatible 1st dimension.");
        if (mxGetDimensions(prhs[4])[1] != dm[1])
            mexErrMsgTxt("Incompatible 2nd dimension.");
        if (mxGetDimensions(prhs[4])[2] != dm[2])
            mexErrMsgTxt("Incompatible 3rd dimension.");
        jd = (float *)mxGetPr(prhs[4]);
    }
    if (mxGetNumberOfElements(prhs[3]) >10)
        mexErrMsgTxt("Fourth argument should contain param1, param2, param3, param4, param5, LMreg, ncycles, nits, nsamps and code.");
    if (mxGetNumberOfElements(prhs[3]) >=1) param[3] = mxGetPr(prhs[3])[0];
    if (mxGetNumberOfElements(prhs[3]) >=2) param[4] = mxGetPr(prhs[3])[1];
    if (mxGetNumberOfElements(prhs[3]) >=3) param[5] = mxGetPr(prhs[3])[2];
    if (mxGetNumberOfElements(prhs[3]) >=4) param[6] = mxGetPr(prhs[3])[3];
    if (mxGetNumberOfElements(prhs[3]) >=5) param[7] = mxGetPr(prhs[3])[4];

    if (mxGetNumberOfElements(prhs[3]) >=6) lmreg0 = mxGetPr(prhs[3])[5];
    if (mxGetNumberOfElements(prhs[3]) >=7) cycles = mxGetPr(prhs[3])[6];
    if (mxGetNumberOfElements(prhs[3]) >=8) its    = mxGetPr(prhs[3])[7];
    if (mxGetNumberOfElements(prhs[3]) >=9) k      = mxGetPr(prhs[3])[8];
    if (mxGetNumberOfElements(prhs[3]) >=10) code  = mxGetPr(prhs[3])[9];

    plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2,nll, mxDOUBLE_CLASS, mxREAL);

    v       = (float *)mxGetPr(prhs[0]);
    g       = (float *)mxGetPr(prhs[1]);
    f       = (float *)mxGetPr(prhs[2]);
    ov      = (float *)mxGetPr(plhs[0]);
    ll      = (double*)mxGetPr(plhs[1]);

    scratch = (float *)mxCalloc(iteration_scratchsize((mwSize *)dm, code,k),sizeof(float));

    dm[3] = 1;
    if (mxGetNumberOfDimensions(prhs[1])>=4)
        dm[3] = mxGetDimensions(prhs[1])[3];

    /* set_bound(0); */
    iteration(dm, k, v, g, f, jd, param, lmreg0, cycles, its, code,
              ov, ll, scratch);
    mxFree((void *)scratch);
}

void exp_mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[])
{
    int k=6;
    mwSize nd;
    const mwSize *dm;
    float *v, *t, *t1;
    double sc = 1.0;
    int flg = 0;

    if (((nrhs != 1) && (nrhs != 2)) || (nlhs>2)) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension must be 3.");

    if (nrhs>1)
    {
        if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[1]) > 3)
            mexErrMsgTxt("Params must contain one to three elements");
        if (mxGetNumberOfElements(prhs[1]) >= 1) k   =   (int)(mxGetPr(prhs[1])[0]);
        if (mxGetNumberOfElements(prhs[1]) >= 2) sc  = (float)(mxGetPr(prhs[1])[1]);
        if (mxGetNumberOfElements(prhs[1]) >= 3) flg =   (int)(mxGetPr(prhs[1])[2]);
    }

    v       = (float *)mxGetPr(prhs[0]);

    plhs[0] = mxCreateNumericArray(nd,dm, mxSINGLE_CLASS, mxREAL);
    t       = (float *)mxGetPr(plhs[0]);
    t1      = mxCalloc(dm[0]*dm[1]*dm[2]*3,sizeof(float));

    /* set_bound(0); */

    if (nlhs < 2)
    {
        expdef((mwSize *)dm, k, sc, v, t, t1, (float *)0, (float *)0);
    }
    else
    {
        float *J, *J1;
        mwSize dmj[5];
        dmj[0]  = dm[0];
        dmj[1]  = dm[1];
        dmj[2]  = dm[2];
        if (flg==0)
        {
            dmj[3]  = 3;
            dmj[4]  = 3;
            plhs[1] = mxCreateNumericArray(5,dmj, mxSINGLE_CLASS, mxREAL);
            J       = (float *)mxGetPr(plhs[1]);
            J1      = mxCalloc(dm[0]*dm[1]*dm[2]*3*3,sizeof(float));
            expdef((mwSize *)dm, k, sc, v, t, t1, J, J1);
        }
        else
        {
            plhs[1] = mxCreateNumericArray(3,dmj, mxSINGLE_CLASS, mxREAL);
            J       = (float *)mxGetPr(plhs[1]);
            J1      = mxCalloc(dm[0]*dm[1]*dm[2],sizeof(float));
            expdefdet((mwSize *)dm, k, sc, v, t, t1, J, J1);
        }
        mxFree((void *)J1);
    }
    unwrap((mwSize *)dm, t);
    mxFree((void *)t1);
}

