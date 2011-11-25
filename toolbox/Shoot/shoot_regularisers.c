/* $Id: shoot_regularisers.c 4573 2011-11-25 23:01:01Z john $ */
/* (c) John Ashburner (2011) */

#include<mex.h>
#include<math.h>
extern double log(double x);

#include "shoot_boundary.h"

#ifdef VERBOSE
static double sumsq_le(mwSize dm[], float a[], float b[], double s[], float u[])
{
    double ss = 0.0;
    mwSignedIndex k;
    double mu = s[3], lam = s[4], id = s[5];
    double wx, wy, wz, wxy, wxz, wyx, wyz, wzx, wzy, w1, w2;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    wx  = 2*mu*(v1+v2)/v0+(4*mu+2*lam) + id/v0;
    wy  = 2*mu*(v0+v2)/v1+(4*mu+2*lam) + id/v1;
    wz  = 2*mu*(v0+v1)/v2+(4*mu+2*lam) + id/v2;
    w1  = -(2*mu+lam);
    wxy = -mu*v0/v1;
    wxz = -mu*v0/v2;
    wyx = -mu*v1/v0;
    wyz = -mu*v1/v2;
    wzx = -mu*v2/v0;
    wzy = -mu*v2/v1;
    w2  = 0.25*(lam+mu);

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km1,kp1;
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *paxy, *payy, *paxz, *payz, *pazz;
            mwSignedIndex i, jm1,jp1;

            pux  = u+dm[0]*(j+dm[1]*k);
            puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
            puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
            pbx  = b+dm[0]*(j+dm[1]*k);
            pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
            pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

            if (a)
            {
                paxx = a+dm[0]*(j+dm[1]*k);
                payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
            }

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im1,ip1;
                float *px = pux+i, *py = puy+i, *pz = puz+i;
                double abx, aby, abz, tmp;

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                if (a)
                {
                    abx = paxx[i]*px[0] + paxy[i]*py[0] + paxz[i]*pz[0];
                    aby = payy[i]*py[0] + paxy[i]*px[0] + payz[i]*pz[0];
                    abz = pazz[i]*pz[0] + paxz[i]*px[0] + payz[i]*py[0];
                }
                else
                    abx = aby = abz = 0.0;

                tmp =  wx * px[0]
                     + w1 *(px[ip1] + px[im1])
                     + wyx*(px[jp1] + px[jm1])
                     + wzx*(px[kp1] + px[km1])
                     + w2 *(py[ip1+jm1] - py[ip1+jp1] - py[im1+jm1] + py[im1+jp1] + pz[ip1+km1] - pz[ip1+kp1] - pz[im1+km1] + pz[im1+kp1])
                     - pbx[i] + abx;
                ss += tmp*tmp;

                tmp =  wy * py[0]
                     + wxy*(py[ip1] + py[im1])
                     + w1 *(py[jp1] + py[jm1])
                     + wzy*(py[kp1] + py[km1])
                     + w2*(px[jp1+im1] - px[jp1+ip1] - px[jm1+im1] + px[jm1+ip1] + pz[jp1+km1] - pz[jp1+kp1] - pz[jm1+km1] + pz[jm1+kp1])
                     - pby[i] + aby;
                ss += tmp*tmp;

                tmp =  wz * pz[0]
                     + wxz*(pz[ip1] + pz[im1])
                     + wyz*(pz[jp1] + pz[jm1])
                     + w1 *(pz[kp1] + pz[km1])
                     + w2 *(px[kp1+im1] - px[kp1+ip1] - px[km1+im1] + px[km1+ip1] + py[kp1+jm1] - py[kp1+jp1] - py[km1+jm1] + py[km1+jp1])
                     - pbz[i] + abz;
                ss += tmp*tmp;
            }
        }
    }
    return(ss);
}
#endif

static void Atimesp1(mwSize dm[], float A[], float p[], float Ap[])
{
    mwSignedIndex i, m = dm[0]*dm[1]*dm[2];
    float *pa11 = A ,     *pa22 = A +m,   *pa33 = A +2*m,
          *pa12 = A +3*m, *pa13 = A +4*m, *pa23 = A +5*m;
    float *pap1 = Ap,     *pap2 = Ap+m,   *pap3 = Ap+2*m;
    float *pp1  = p,      *pp2  = p +m,   *pp3  = p +2*m;

    if (A==0) return;

    for(i=0; i<m; i++)
    {
        pap1[i] += pa11[i]*pp1[i] + pa12[i]*pp2[i] + pa13[i]*pp3[i];
        pap2[i] += pa12[i]*pp1[i] + pa22[i]*pp2[i] + pa23[i]*pp3[i];
        pap3[i] += pa13[i]*pp1[i] + pa23[i]*pp2[i] + pa33[i]*pp3[i];
    }
}


void vel2mom_le(mwSize dm[], float f[], double s[], float g[])
{
    mwSignedIndex k;
    double mu = s[3], lam = s[4], id = s[5];
    double wx, wy, wz, wxy, wxz, wyx, wyz, wzx, wzy, w1, w2;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    wx  = 2*mu*(v1+v2)/v0+(4*mu+2*lam) + id/v0;
    wy  = 2*mu*(v0+v2)/v1+(4*mu+2*lam) + id/v1;
    wz  = 2*mu*(v0+v1)/v2+(4*mu+2*lam) + id/v2;
    w1  = -(2*mu+lam);
    wxy = -mu*v0/v1;
    wxz = -mu*v0/v2;
    wyx = -mu*v1/v0;
    wyz = -mu*v1/v2;
    wzx = -mu*v2/v0;
    wzy = -mu*v2/v1;
    w2  = 0.25*(lam+mu);

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km1,kp1;
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex i, jm1,jp1;
            float *pgx, *pgy, *pgz, *pfx, *pfy, *pfz;

            pgx = g+dm[0]*(j+dm[1]*k);
            pgy = g+dm[0]*(j+dm[1]*(k+dm[2]));
            pgz = g+dm[0]*(j+dm[1]*(k+dm[2]*2));

            pfx = f+dm[0]*(j+dm[1]*k);
            pfy = f+dm[0]*(j+dm[1]*(k+dm[2]));
            pfz = f+dm[0]*(j+dm[1]*(k+dm[2]*2));

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im1,ip1;
                float *px = &pfx[i], *py = &pfy[i], *pz = &pfz[i];

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                pgx[i] = wx *px[0]
                       + w1 *(px[ip1] + px[im1])
                       + wyx*(px[jp1] + px[jm1])
                       + wzx*(px[kp1] + px[km1])
                       + w2 *(py[ip1+jm1] - py[ip1+jp1] - py[im1+jm1] + py[im1+jp1] + pz[ip1+km1] - pz[ip1+kp1] - pz[im1+km1] + pz[im1+kp1]);

                pgy[i] = wy *py[0]
                       + wxy*(py[ip1] + py[im1])
                       + w1 *(py[jp1] + py[jm1])
                       + wzy*(py[kp1] + py[km1])
                       + w2 *(px[jp1+im1] - px[jp1+ip1] - px[jm1+im1] + px[jm1+ip1] + pz[jp1+km1] - pz[jp1+kp1] - pz[jm1+km1] + pz[jm1+kp1]);

                pgz[i] = wz *pz[0]
                       + wxz*(pz[ip1] + pz[im1])
                       + wyz*(pz[jp1] + pz[jm1])
                       + w1 *(pz[kp1] + pz[km1])
                       + w2 *(px[kp1+im1] - px[kp1+ip1] - px[km1+im1] + px[km1+ip1] + py[kp1+jm1] - py[kp1+jp1] - py[km1+jm1] + py[km1+jp1]);
            }
        }
    }
}

void relax_le(mwSize dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double mu = s[3], lam = s[4], id = s[5];
    double wx, wy, wz, wxy, wxz, wyx, wyz, wzx, wzy, w1, w2;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    wx  = 2*mu*(v1+v2)/v0+(4*mu+2*lam) + id/v0;
    wy  = 2*mu*(v0+v2)/v1+(4*mu+2*lam) + id/v1;
    wz  = 2*mu*(v0+v1)/v2+(4*mu+2*lam) + id/v2;
    w1  = -(2*mu+lam);
    wxy = -mu*v0/v1;
    wxz = -mu*v0/v2;
    wyx = -mu*v1/v0;
    wyz = -mu*v1/v2;
    wzx = -mu*v2/v0;
    wzy = -mu*v2/v1;
    w2  = 0.25*(lam+mu);

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("%dx%dx%d: ", dm[0],dm[1],dm[2]);
#   endif

    for(it=0; it<8*nit; it++)
    {
        mwSignedIndex k;
        for(k=it&1; k<dm[2]; k+=2)
        {
            mwSignedIndex j, km1, kp1;
            km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

            for(j=(it>>1)&1; j<dm[1]; j+=2)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
                mwSignedIndex i, jm1,jp1;

                pux  = u+dm[0]*(j+dm[1]* k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]* k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

                if (a)
                {
                    paxx = a+dm[0]*(j+dm[1]* k);
                    payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                    pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                    paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                    paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                    payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
                }

                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

                for(i=(it>>2)&1; i<dm[0]; i+=2)
                {
                    mwSignedIndex im1,ip1;
                    double sux, suy, suz;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;

                    sux = pbx[i] - ( w1 *(px[ip1] + px[im1])
                                   + wyx*(px[jp1] + px[jm1])
                                   + wzx*(px[kp1] + px[km1])
                                   + w2 *(py[ip1+jm1] - py[ip1+jp1] - py[im1+jm1] + py[im1+jp1] + pz[ip1+km1] - pz[ip1+kp1] - pz[im1+km1] + pz[im1+kp1]));

                    suy = pby[i] - ( wxy*(py[ip1] + py[im1])
                                   + w1 *(py[jp1] + py[jm1])
                                   + wzy*(py[kp1] + py[km1])
                                   + w2 *(px[jp1+im1] - px[jp1+ip1] - px[jm1+im1] + px[jm1+ip1] + pz[jp1+km1] - pz[jp1+kp1] - pz[jm1+km1] + pz[jm1+kp1]));

                    suz = pbz[i] - ( wxz*(pz[ip1] + pz[im1])
                                   + wyz*(pz[jp1] + pz[jm1])
                                   + w1 *(pz[kp1] + pz[km1])
                                   + w2 *(px[kp1+im1] - px[kp1+ip1] - px[km1+im1] + px[km1+ip1] + py[kp1+jm1] - py[kp1+jp1] - py[km1+jm1] + py[km1+jp1]));
                    if (a)
                    {
                        double axx, ayy, azz, axy, axz, ayz, idt;

                        axx  = paxx[i] + wx;
                        ayy  = payy[i] + wy;
                        azz  = pazz[i] + wz;
                        axy  = paxy[i];
                        axz  = paxz[i];
                        ayz  = payz[i];
                        idt  = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);

                        *px = idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                        *py = idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                        *pz = idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                    }
                    else
                    {
                        *px = sux/wx;
                        *py = suy/wy;
                        *pz = suz/wz;
                    }
                }
            }
        } 
#       ifdef VERBOSE
            if ((it%8)==7) printf(" %g", sumsq_le(dm, a, b, s, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}


void Atimesp_le(mwSize dm[], float A[], double param[], float p[], float Ap[])
{
    vel2mom_le(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

#ifdef VERBOSE
static double sumsq_me(mwSize dm[], float a[], float b[], double s[], float u[])
{
    double w000, w001, w010, w100;
    double ss = 0.0;
    mwSignedIndex i, j, k;

    w000 = s[3]*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + s[5];
    w001 = s[3]*(-s[2]*s[2]);
    w010 = s[3]*(-s[1]*s[1]);
    w100 = s[3]*(-s[0]*s[0]);

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex km1,kp1;
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *paxy, *payy, *paxz, *payz, *pazz;
            mwSignedIndex jm1,jp1,im1,ip1;

            pux  = u+dm[0]*(j+dm[1]*k);
            puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
            puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
            pbx  = b+dm[0]*(j+dm[1]*k);
            pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
            pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

            if (a)
            {
                paxx = a+dm[0]*(j+dm[1]*k);
                payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
            }

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                float *px = &pux[i], *py = &puy[i], *pz = &puz[i];
                double abx, aby, abz, tmp;

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                if (a)
                {
                    abx = paxx[i]*px[0] + paxy[i]*py[0] + paxz[i]*pz[0];
                    aby = paxy[i]*px[0] + payy[i]*py[0] + payz[i]*pz[0];
                    abz = paxz[i]*px[0] + payz[i]*py[0] + pazz[i]*pz[0];
                }
                else
                {
                    abx = aby = abz = 0.0;
                }
                tmp =  (w000* px[0] 
                      + w001*(px[km1] + px[kp1])
                      + w010*(px[jm1] + px[jp1])
                      + w100*(px[im1] + px[ip1]))/(s[0]*s[0])
                      - pbx[i] + abx;
                ss += tmp*tmp;

                tmp =  (w000* py[0]
                      + w001*(py[km1] + py[kp1])
                      + w010*(py[jm1] + py[jp1])
                      + w100*(py[im1] + py[ip1]))/(s[1]*s[1])
                      - pby[i] + aby;
                ss += tmp*tmp;

                tmp =  (w000* pz[0]
                      + w001*(pz[km1] + pz[kp1])
                      + w010*(pz[jm1] + pz[jp1])
                      + w100*(pz[im1] + pz[ip1]))/(s[2]*s[2])
                      - pbz[i] + abz;
                ss += tmp*tmp;
            }
        }
    }
    return(ss);
}
#endif

void vel2mom_me(mwSize dm[], float f[], double s[], float g[])
{
    mwSignedIndex i, j, k, km1,kp1, jm1,jp1, im1,ip1;
    float *pgx, *pgy, *pgz, *pfx, *pfy, *pfz;
    double w000,w001,w010,w100;

    w000 = s[3]*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + s[5];
    w001 = s[3]*(-s[2]*s[2]);
    w010 = s[3]*(-s[1]*s[1]);
    w100 = s[3]*(-s[0]*s[0]);

    for(k=0; k<dm[2]; k++)
    {
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            pgx = g+dm[0]*(j+dm[1]*k);
            pgy = g+dm[0]*(j+dm[1]*(k+dm[2]));
            pgz = g+dm[0]*(j+dm[1]*(k+dm[2]*2));

            pfx = f+dm[0]*(j+dm[1]*k);
            pfy = f+dm[0]*(j+dm[1]*(k+dm[2]));
            pfz = f+dm[0]*(j+dm[1]*(k+dm[2]*2));

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                float *px = &pfx[i], *py = &pfy[i], *pz = &pfz[i];

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                pgx[i] = (w000*px[0] + w001*(px[km1] + px[kp1]) + w010*(px[jm1] + px[jp1]) + w100*(px[im1] + px[ip1]))/(s[0]*s[0]);
                pgy[i] = (w000*py[0] + w001*(py[km1] + py[kp1]) + w010*(py[jm1] + py[jp1]) + w100*(py[im1] + py[ip1]))/(s[1]*s[1]);
                pgz[i] = (w000*pz[0] + w001*(pz[km1] + pz[kp1]) + w010*(pz[jm1] + pz[jp1]) + w100*(pz[im1] + pz[ip1]))/(s[2]*s[2]);
            }
        }
    }
}

void relax_me(mwSize dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double w000,w001,w010,w100;

    w000 = s[3]*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + s[5];
    w001 = s[3]*(-s[2]*s[2]);
    w010 = s[3]*(-s[1]*s[1]);
    w100 = s[3]*(-s[0]*s[0]);

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("%dx%dx%d: ", dm[0],dm[1],dm[2]);
#   endif

    for(it=0; it<2*nit; it++)
    {
        mwSignedIndex k, kstart;
        mwSignedIndex j, jstart;
        mwSignedIndex i, istart;

        kstart = it%2;
        for(k=0; k<dm[2]; k++)
        {
            mwSignedIndex km1, kp1;
            km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

            jstart = (kstart == (k%2));
            for(j=0; j<dm[1]; j++)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *paxy, *payy, *paxz, *payz, *pazz;
                mwSignedIndex jm1,jp1, im1,ip1;

                pux  = u+dm[0]*(j+dm[1]*k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]*k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

                if (a)
                {
                    paxx = a+dm[0]*(j+dm[1]*k);
                    payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                    pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                    paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                    paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                    payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
                }

                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

                istart = (jstart == (j%2));

                for(i=istart; i<dm[0]; i+=2)
                {
                    double sux, suy, suz;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;

                    sux = pbx[i]-(w001*(px[km1] + px[kp1]) + w010*(px[jm1] + px[jp1]) + w100*(px[im1] + px[ip1]))/(s[0]*s[0]);
                    suy = pby[i]-(w001*(py[km1] + py[kp1]) + w010*(py[jm1] + py[jp1]) + w100*(py[im1] + py[ip1]))/(s[1]*s[1]);
                    suz = pbz[i]-(w001*(pz[km1] + pz[kp1]) + w010*(pz[jm1] + pz[jp1]) + w100*(pz[im1] + pz[ip1]))/(s[2]*s[2]);

                    if (a)
                    {
                        double axx, ayy, azz, axy, axz, ayz, idt;
                        /*
                           syms axx ayy azz axy axz ayz sux suy suz
                           A = [axx axy axz; axy ayy ayz; axz ayz azz];
                           su = [sux ; suy; suz]
                           simplify(inv(A)*su)
                        */
                        axx = paxx[i] + w000/(s[0]*s[0]);
                        ayy = payy[i] + w000/(s[1]*s[1]);
                        azz = pazz[i] + w000/(s[2]*s[2]);
                        axy = paxy[i];
                        axz = paxz[i];
                        ayz = payz[i];
                        idt = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);
                        *px = idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                        *py = idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                        *pz = idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                    }
                    else
                    {
                        *px = (s[0]*s[0])*sux/w000;
                        *py = (s[1]*s[1])*suy/w000;
                        *pz = (s[2]*s[2])*suz/w000;
                    }
                }
            }
        }
#   ifdef VERBOSE
        if ((it%2)==1) printf(" %g", sumsq_me(dm, a, b, s, u));
#   endif
    }
#ifdef VERBOSE
    printf("\n");
#endif

}

void Atimesp_me(mwSize dm[], float A[], double param[], float p[], float Ap[])
{
    vel2mom_me(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

#ifdef VERBOSE
static double sumsq_be(mwSize dm[], float a[], float b[], double s[], float u[])
{
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double ss = 0.0;
    mwSignedIndex k;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    w000 = s[3]*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +s[4]*2*(v0+v1+v2) + s[5];
    w100 = s[3]*(-4*v0*(v0+v1+v2)) -s[4]*v0;
    w010 = s[3]*(-4*v1*(v0+v1+v2)) -s[4]*v1;
    w001 = s[3]*(-4*v2*(v0+v1+v2)) -s[4]*v2;
    w200 = s[3]*v0*v0;
    w020 = s[3]*v1*v1;
    w002 = s[3]*v2*v2;
    w110 = s[3]*2*v0*v1;
    w101 = s[3]*2*v0*v2;
    w011 = s[3]*2*v1*v2;

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km2,km1,kp1,kp2;
        km2 = (BOUND(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (BOUND(k+2,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
            mwSignedIndex i, jm2,jm1,jp1,jp2;

            pux  = u+dm[0]*(j+dm[1]*k);
            puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
            puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
            pbx  = b+dm[0]*(j+dm[1]*k);
            pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
            pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

            if (a)
            {
                paxx = a+dm[0]*(j+dm[1]*k);
                payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
            }

            jm2 = (BOUND(j-2,dm[1])-j)*dm[0];
            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];
            jp2 = (BOUND(j+2,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im2,im1,ip1,ip2;
                float *px = pux+i, *py = puy+i, *pz = puz+i;
                double tmp, abx, aby, abz;

                im2 = BOUND(i-2,dm[0])-i;
                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;
                ip2 = BOUND(i+2,dm[0])-i;

                if (a)
                {
                    abx = paxx[i]*px[0] + paxy[i]*py[0] + paxz[i]*pz[0];
                    aby = paxy[i]*px[0] + payy[i]*py[0] + payz[i]*pz[0];
                    abz = paxz[i]*px[0] + payz[i]*py[0] + pazz[i]*pz[0];
                }
                else
                {
                    abx = aby = abz = 0.0;
                }
                tmp = abx
                     +(w000* px[0]
                     + w010*(px[    jm1    ] + px[    jp1    ])
                     + w020*(px[    jm2    ] + px[    jp2    ])
                     + w100*(px[im1        ] + px[ip1        ])
                     + w110*(px[im1+jm1    ] + px[ip1+jm1    ] + px[im1+jp1    ] + px[ip1+jp1    ])
                     + w200*(px[im2        ] + px[ip2        ])
                     + w001*(px[        km1] + px[        kp1])
                     + w101*(px[im1    +km1] + px[ip1    +km1] + px[im1    +kp1] + px[ip1    +kp1])
                     + w011*(px[    jm1+km1] + px[    jp1+km1] + px[    jm1+kp1] + px[    jp1+kp1])
                     + w002*(px[        km2] + px[        kp2]))/v0
                     - pbx[i];
                ss += tmp*tmp;

                tmp = aby
                     +(w000* py[0]
                     + w010*(py[    jm1    ] + py[    jp1    ])
                     + w020*(py[    jm2    ] + py[    jp2    ])
                     + w100*(py[im1        ] + py[ip1        ])
                     + w110*(py[im1+jm1    ] + py[ip1+jm1    ] + py[im1+jp1    ] + py[ip1+jp1    ])
                     + w200*(py[im2        ] + py[ip2        ])
                     + w001*(py[        km1] + py[        kp1])
                     + w101*(py[im1    +km1] + py[ip1    +km1] + py[im1    +kp1] + py[ip1    +kp1])
                     + w011*(py[    jm1+km1] + py[    jp1+km1] + py[    jm1+kp1] + py[    jp1+kp1])
                     + w002*(py[        km2] + py[        kp2]))/v1
                     - pby[i];
                ss += tmp*tmp;

                tmp = abz
                     +(w000* pz[0]
                     + w010*(pz[    jm1    ] + pz[    jp1    ])
                     + w020*(pz[    jm2    ] + pz[    jp2    ])
                     + w100*(pz[im1        ] + pz[ip1        ])
                     + w110*(pz[im1+jm1    ] + pz[ip1+jm1    ] + pz[im1+jp1    ] + pz[ip1+jp1    ])
                     + w200*(pz[im2        ] + pz[ip2        ])
                     + w001*(pz[        km1] + pz[        kp1])
                     + w101*(pz[im1    +km1] + pz[ip1    +km1] + pz[im1    +kp1] + pz[ip1    +kp1])
                     + w011*(pz[    jm1+km1] + pz[    jp1+km1] + pz[    jm1+kp1] + pz[    jp1+kp1])
                     + w002*(pz[        km2] + pz[        kp2]))/v2
                     - pbz[i];
                ss += tmp*tmp;
            }
        }
    }
    return(ss);
}
#endif

void vel2mom_be(mwSize dm[], float f[], double s[], float g[])
{
    mwSignedIndex k;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    w000 = s[3]*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +s[4]*2*(v0+v1+v2) + s[5];
    w100 = s[3]*(-4*v0*(v0+v1+v2)) -s[4]*v0;
    w010 = s[3]*(-4*v1*(v0+v1+v2)) -s[4]*v1;
    w001 = s[3]*(-4*v2*(v0+v1+v2)) -s[4]*v2;
    w200 = s[3]*v0*v0;
    w020 = s[3]*v1*v1;
    w002 = s[3]*v2*v2;
    w110 = s[3]*2*v0*v1;
    w101 = s[3]*2*v0*v2;
    w011 = s[3]*2*v1*v2;

    /*
        syms s1 s2 s3
        syms l1 l2 l3
        zz = sym(zeros(3,3));
        K1 = cat(3,zz,[0 -s1*s1 0; 0 2*s1*s1 0; 0 -s1*s1 0],zz);
        K2 = cat(3,zz,[0 0 0; -s2*s2 2*s2*s2 -s2*s2; 0 0 0],zz);
        K3 = sym(zeros(3,3,3));
        K3(2,2,1) = -s3*s3;
        K3(2,2,2) = 2*s3*s3;
        K3(2,2,3) = -s3*s3;

        K  = K1+K2+K3;
        K1 = K*l1; K1(2,2,2) = K1(2,2,2)+l2;
        K2 = K*l1; K2(2,2,2) = K2(2,2,2)+l3;

        % L  = convn(K,K)
        L  = sym(zeros(5,5,5));
        for i=1:3,
            for j=1:3,
                for k=1:3,
                    L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) = L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) + K1(i,j,k)*K2;
                end;
            end;
        end;
        disp(simplify(L(3:end,3:end,3:end)))
    */

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km2,km1,kp1,kp2;
        km2 = (BOUND(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (BOUND(k+2,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex i, jm2,jm1,jp1,jp2;
            float *pgx, *pgy, *pgz, *pfx, *pfy, *pfz;

            pgx = g+dm[0]*(j+dm[1]*k);
            pgy = g+dm[0]*(j+dm[1]*(k+dm[2]));
            pgz = g+dm[0]*(j+dm[1]*(k+dm[2]*2));

            pfx = f+dm[0]*(j+dm[1]*k);
            pfy = f+dm[0]*(j+dm[1]*(k+dm[2]));
            pfz = f+dm[0]*(j+dm[1]*(k+dm[2]*2));

            jm2 = (BOUND(j-2,dm[1])-j)*dm[0];
            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];
            jp2 = (BOUND(j+2,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im2,im1,ip1,ip2;
                float *px = &pfx[i], *py = &pfy[i], *pz = &pfz[i];

                im2 = BOUND(i-2,dm[0])-i;
                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;
                ip2 = BOUND(i+2,dm[0])-i;

                pgx[i] =(w000* px[0]
                       + w010*(px[    jm1    ] + px[    jp1    ])
                       + w020*(px[    jm2    ] + px[    jp2    ])
                       + w100*(px[im1        ] + px[ip1        ])
                       + w110*(px[im1+jm1    ] + px[ip1+jm1    ] + px[im1+jp1    ] + px[ip1+jp1    ])
                       + w200*(px[im2        ] + px[ip2        ])
                       + w001*(px[        km1] + px[        kp1])
                       + w101*(px[im1    +km1] + px[ip1    +km1] + px[im1    +kp1] + px[ip1    +kp1])
                       + w011*(px[    jm1+km1] + px[    jp1+km1] + px[    jm1+kp1] + px[    jp1+kp1])
                       + w002*(px[        km2] + px[        kp2]))/v0;
                pgy[i] =(w000* py[0]
                       + w010*(py[    jm1    ] + py[    jp1    ])
                       + w020*(py[    jm2    ] + py[    jp2    ])
                       + w100*(py[im1        ] + py[ip1        ])
                       + w110*(py[im1+jm1    ] + py[ip1+jm1    ] + py[im1+jp1    ] + py[ip1+jp1    ])
                       + w200*(py[im2        ] + py[ip2        ])
                       + w001*(py[        km1] + py[        kp1])
                       + w101*(py[im1    +km1] + py[ip1    +km1] + py[im1    +kp1] + py[ip1    +kp1])
                       + w011*(py[    jm1+km1] + py[    jp1+km1] + py[    jm1+kp1] + py[    jp1+kp1])
                       + w002*(py[        km2] + py[        kp2]))/v1;
                pgz[i] =(w000* pz[0]
                       + w010*(pz[    jm1    ] + pz[    jp1    ])
                       + w020*(pz[    jm2    ] + pz[    jp2    ])
                       + w100*(pz[im1        ] + pz[ip1        ])
                       + w110*(pz[im1+jm1    ] + pz[ip1+jm1    ] + pz[im1+jp1    ] + pz[ip1+jp1    ])
                       + w200*(pz[im2        ] + pz[ip2        ])
                       + w001*(pz[        km1] + pz[        kp1])
                       + w101*(pz[im1    +km1] + pz[ip1    +km1] + pz[im1    +kp1] + pz[ip1    +kp1])
                       + w011*(pz[    jm1+km1] + pz[    jp1+km1] + pz[    jm1+kp1] + pz[    jp1+kp1])
                       + w002*(pz[        km2] + pz[        kp2]))/v2;
            }
        }
    }
}

void relax_be(mwSize dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    w000 = s[3]*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +s[4]*2*(v0+v1+v2) + s[5];
    w100 = s[3]*(-4*v0*(v0+v1+v2)) -s[4]*v0;
    w010 = s[3]*(-4*v1*(v0+v1+v2)) -s[4]*v1;
    w001 = s[3]*(-4*v2*(v0+v1+v2)) -s[4]*v2;
    w200 = s[3]*v0*v0;
    w020 = s[3]*v1*v1;
    w002 = s[3]*v2*v2;
    w110 = s[3]*2*v0*v1;
    w101 = s[3]*2*v0*v2;
    w011 = s[3]*2*v1*v2;

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("%dx%dx%d: ", dm[0],dm[1],dm[2]);
#   endif

    for(it=0; it<27*nit; it++)
    {
        mwSignedIndex i, j, k;
        for(k=(it/9)%3; k<dm[2]; k+=3)
        {
            mwSignedIndex km2, km1, kp1, kp2;
            km2 = (BOUND(k-2,dm[2])-k)*dm[0]*dm[1];
            km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];
            kp2 = (BOUND(k+2,dm[2])-k)*dm[0]*dm[1];

            for(j=(it/3)%3; j<dm[1]; j+=3)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
                mwSignedIndex jm2,jm1,jp1,jp2;

                pux  = u+dm[0]*(j+dm[1]* k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]* k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

                if (a)
                {
                    paxx = a+dm[0]*(j+dm[1]* k);
                    payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                    pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                    paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                    paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                    payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
                }

                jm2 = (BOUND(j-2,dm[1])-j)*dm[0];
                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];
                jp2 = (BOUND(j+2,dm[1])-j)*dm[0];

                for(i=it%3; i<dm[0]; i+=3)
                {
                    mwSignedIndex im2,im1,ip1,ip2;
                    double sux, suy, suz;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im2 = BOUND(i-2,dm[0])-i;
                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;
                    ip2 = BOUND(i+2,dm[0])-i;

                    sux = pbx[i] - (w010*(px[    jm1    ] + px[    jp1    ])
                                  + w020*(px[    jm2    ] + px[    jp2    ])
                                  + w100*(px[im1        ] + px[ip1        ])
                                  + w110*(px[im1+jm1    ] + px[ip1+jm1    ] + px[im1+jp1    ] + px[ip1+jp1    ])
                                  + w200*(px[im2        ] + px[ip2        ])
                                  + w001*(px[        km1] + px[        kp1])
                                  + w101*(px[im1    +km1] + px[ip1    +km1] + px[im1    +kp1] + px[ip1    +kp1])
                                  + w011*(px[    jm1+km1] + px[    jp1+km1] + px[    jm1+kp1] + px[    jp1+kp1])
                                  + w002*(px[        km2] + px[        kp2]))/v0;

                    suy = pby[i] - (w010*(py[    jm1    ] + py[    jp1    ])
                                  + w020*(py[    jm2    ] + py[    jp2    ])
                                  + w100*(py[im1        ] + py[ip1        ])
                                  + w110*(py[im1+jm1    ] + py[ip1+jm1    ] + py[im1+jp1    ] + py[ip1+jp1    ])
                                  + w200*(py[im2        ] + py[ip2        ])
                                  + w001*(py[        km1] + py[        kp1])
                                  + w101*(py[im1    +km1] + py[ip1    +km1] + py[im1    +kp1] + py[ip1    +kp1])
                                  + w011*(py[    jm1+km1] + py[    jp1+km1] + py[    jm1+kp1] + py[    jp1+kp1])
                                  + w002*(py[        km2] + py[        kp2]))/v1;

                    suz = pbz[i] - (w010*(pz[    jm1    ] + pz[    jp1    ])
                                  + w020*(pz[    jm2    ] + pz[    jp2    ])
                                  + w100*(pz[im1        ] + pz[ip1        ])
                                  + w110*(pz[im1+jm1    ] + pz[ip1+jm1    ] + pz[im1+jp1    ] + pz[ip1+jp1    ])
                                  + w200*(pz[im2        ] + pz[ip2        ])
                                  + w001*(pz[        km1] + pz[        kp1])
                                  + w101*(pz[im1    +km1] + pz[ip1    +km1] + pz[im1    +kp1] + pz[ip1    +kp1])
                                  + w011*(pz[    jm1+km1] + pz[    jp1+km1] + pz[    jm1+kp1] + pz[    jp1+kp1])
                                  + w002*(pz[        km2] + pz[        kp2]))/v2;

                    if (a)
                    {
                        double axx, ayy, azz, axy, axz, ayz, idt;
/*
                        sux -= (paxx[i]*px[0] + paxy[i]*py[0] + paxz[i]*pz[0]);
                        suy -= (paxy[i]*px[0] + payy[i]*py[0] + payz[i]*pz[0]);
                        suz -= (paxz[i]*px[0] + payz[i]*py[0] + pazz[i]*pz[0]);
*/
                        axx  = paxx[i] + w000/v0;
                        ayy  = payy[i] + w000/v1;
                        azz  = pazz[i] + w000/v2;
                        axy  = paxy[i];
                        axz  = paxz[i];
                        ayz  = payz[i];
                        idt  = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);
                        *px = idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                        *py = idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                        *pz = idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                    }
                    else
                    {
                        *px = v0*sux/w000;
                        *py = v1*suy/w000;
                        *pz = v2*suz/w000;
                    }
                }
            }
        }
#       ifdef VERBOSE
        if ((it%27) == 26)
            printf(" %g", sumsq_be(dm, a, b, s, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}


void Atimesp_be(mwSize dm[], float A[], double param[], float p[], float Ap[])
{
    vel2mom_be(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

