/* $Id: optimizer3d.c 4016 2010-07-26 13:12:40Z john $ */
/* (c) John Ashburner (2007) */

#include<mex.h>
#include<math.h>
extern double log(double x);

#include "optimizer3d.h"

#ifdef NEUMANN
    /* Neumann boundary condition */
    static int neumann(int i, int m)
    {
        if (m==1)
            return(0);
        else
        {
            int m2 = m*2;
            i = (i<0) ? (-i-m2*((-i)/m2)-1) : (i-m2*(i/m2));
            if (m<=i)
                return(m2-i-1);
            else
                return(i);
        }
    }
#   define BOUND(i,m) neumann(i,m)
#else
    /* circulant boundary condition */
#   define BOUND(i,m) (((i)>=0) ? (i)%(m) : ((m)+(i)%(m))%m)
#endif

static double sumsq_le_noa(int dm[], float b[], double s[], float u[])
{
    double ss = 0.0;
    int k;
    double mu = s[3], lam = s[4], id = s[5];
    double wx0, wx1, wx2, wy0, wy1, wy2, wz0, wz1, wz2, wxy, wxz, wyz;

    wx0 = 2*mu*(s[1]*s[1]+s[2]*s[2])+(4*mu+2*lam)*s[0]*s[0] + id;
    wy0 = 2*mu*(s[0]*s[0]+s[2]*s[2])+(4*mu+2*lam)*s[1]*s[1] + id;
    wz0 = 2*mu*(s[0]*s[0]+s[1]*s[1])+(4*mu+2*lam)*s[2]*s[2] + id;

    wx1 = -(2*mu+lam)*s[0]*s[0];
    wy1 = -(2*mu+lam)*s[1]*s[1];
    wz1 = -(2*mu+lam)*s[2]*s[2];

    wx2 = -mu*s[0]*s[0];
    wy2 = -mu*s[1]*s[1];
    wz2 = -mu*s[2]*s[2];

    wxy = 0.25*(lam+mu)*s[0]*s[1];
    wxz = 0.25*(lam+mu)*s[0]*s[2];
    wyz = 0.25*(lam+mu)*s[1]*s[2];

    for(k=0; k<dm[2]; k++)
    {
        int j, km1,kp1;
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            float *pux, *puy, *puz, *pbx, *pby, *pbz;
            int i, jm1,jp1;

            pux  = u+dm[0]*(j+dm[1]*k);
            puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
            puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
            pbx  = b+dm[0]*(j+dm[1]*k);
            pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
            pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                int im1,ip1;
                float *px = pux+i, *py = puy+i, *pz = puz+i;
                double tmp;

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                tmp =  wx0*px[0]
                     + wx1*(px[ip1] + px[im1])
                     + wy2*(px[jp1] + px[jm1])
                     + wz2*(px[kp1] + px[km1])
                     + wxy*(py[ip1+jm1] - py[ip1+jp1] - py[im1+jm1] + py[im1+jp1])
                     + wxz*(pz[ip1+km1] - pz[ip1+kp1] - pz[im1+km1] + pz[im1+kp1])
                     - pbx[i];
                ss += tmp*tmp;

                tmp =  wy0*py[0]
                     + wx2*(py[ip1] + py[im1])
                     + wy1*(py[jp1] + py[jm1])
                     + wz2*(py[kp1] + py[km1])
                     + wxy*(px[jp1+im1] - px[jp1+ip1] - px[jm1+im1] + px[jm1+ip1])
                     + wyz*(pz[jp1+km1] - pz[jp1+kp1] - pz[jm1+km1] + pz[jm1+kp1])
                     - pby[i];
                ss += tmp*tmp;

                tmp =  wz0*pz[0]
                     + wx2*(pz[ip1] + pz[im1])
                     + wy2*(pz[jp1] + pz[jm1])
                     + wz1*(pz[kp1] + pz[km1])
                     + wxz*(px[kp1+im1] - px[kp1+ip1] - px[km1+im1] + px[km1+ip1])
                     + wyz*(py[kp1+jm1] - py[kp1+jp1] - py[km1+jm1] + py[km1+jp1])
                     - pbz[i];
                ss += tmp*tmp;

            }
        }
    }
    return(ss);
}

static void relax_le_noa(int dm[], float b[], double s[], int nit, float u[])
{
    int it;
    double regx, regy, regz;
    double mu = s[3], lam = s[4], id = s[5];
    double wx0, wx1, wx2, wy0, wy1, wy2, wz0, wz1, wz2, wxy, wxz, wyz;

    wx0 = 2*mu*(s[1]*s[1]+s[2]*s[2])+(4*mu+2*lam)*s[0]*s[0] + id;
    wy0 = 2*mu*(s[0]*s[0]+s[2]*s[2])+(4*mu+2*lam)*s[1]*s[1] + id;
    wz0 = 2*mu*(s[0]*s[0]+s[1]*s[1])+(4*mu+2*lam)*s[2]*s[2] + id;

    wx1 = -(2*mu+lam)*s[0]*s[0];
    wy1 = -(2*mu+lam)*s[1]*s[1];
    wz1 = -(2*mu+lam)*s[2]*s[2];

    wx2 = -mu*s[0]*s[0];
    wy2 = -mu*s[1]*s[1];
    wz2 = -mu*s[2]*s[2];

    wxy = 0.25*(lam+mu)*s[0]*s[1];
    wxz = 0.25*(lam+mu)*s[0]*s[2];
    wyz = 0.25*(lam+mu)*s[1]*s[2];

    /* For stability in Gauss-Seidel relaxation, the magnitude of the diagonal element must
       exceed the sum of the magnitudes of the off diagonal elements of each column or row
       (see e.g. http://www.mathpages.com/home/kmath175/kmath175.htm).
       This should stabilise the relaxation, providing the second derives are positive definite. */
    regx = (4.0*(wxy+wxz)-2.0*(wx1+wy2+wz2)) - wx0; if (regx<0.0) regx = 0.0; regx += wx0;
    regy = (4.0*(wxy+wyz)-2.0*(wx2+wy1+wz2)) - wy0; if (regy<0.0) regy = 0.0; regy += wy0;
    regz = (4.0*(wxz+wyz)-2.0*(wx2+wy2+wz1)) - wz0; if (regz<0.0) regz = 0.0; regz += wz0;

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("%dx%dx%d (%g,%g,%g): ", dm[0],dm[1],dm[2], regx,regy,regz);
#   endif

    for(it=0; it<8*nit; it++)
    {
        int k;
        /* double ss = 0.0; */
        for(k=it&1; k<dm[2]; k+=2)
        {
            int j, km1, kp1;
            km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

            for(j=(it>>1)&1; j<dm[1]; j+=2)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz;
                int i, jm1,jp1;

                pux  = u+dm[0]*(j+dm[1]* k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]* k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

                for(i=(it>>2)&1; i<dm[0]; i+=2)
                {
                    int im1,ip1;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;

                    *px +=(pbx[i] - (wx0*px[0]
                                   + wx1*(px[ip1] + px[im1])
                                   + wy2*(px[jp1] + px[jm1])
                                   + wz2*(px[kp1] + px[km1])
                                   + wxy*(py[ip1+jm1] - py[ip1+jp1] - py[im1+jm1] + py[im1+jp1])
                                   + wxz*(pz[ip1+km1] - pz[ip1+kp1] - pz[im1+km1] + pz[im1+kp1])))/regx;;

                    *py +=(pby[i] - (wy0*py[0]
                                   + wx2*(py[ip1] + py[im1])
                                   + wy1*(py[jp1] + py[jm1])
                                   + wz2*(py[kp1] + py[km1])
                                   + wxy*(px[jp1+im1] - px[jp1+ip1] - px[jm1+im1] + px[jm1+ip1])
                                   + wyz*(pz[jp1+km1] - pz[jp1+kp1] - pz[jm1+km1] + pz[jm1+kp1])))/regy;

                    *pz +=(pbz[i] - (wz0*pz[0]
                                   + wx2*(pz[ip1] + pz[im1])
                                   + wy2*(pz[jp1] + pz[jm1])
                                   + wz1*(pz[kp1] + pz[km1])
                                   + wxz*(px[kp1+im1] - px[kp1+ip1] - px[km1+im1] + px[km1+ip1])
                                   + wyz*(py[kp1+jm1] - py[kp1+jp1] - py[km1+jm1] + py[km1+jp1])))/regz;
                }
            }
        }

#       ifdef VERBOSE
            printf(" %g", sumsq_le_noa(dm, b, s, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}

static void solve33_noa(float b[], double t, float u[])
{
    u[0] = b[0]/t;
    u[1] = b[1]/t;
    u[2] = b[2]/t;
}

static void Atimesp1(int dm[], float A[], float p[], float Ap[])
{
    int i, m = dm[0]*dm[1]*dm[2];
    float *pa11 = A ,     *pa22 = A +m,   *pa33 = A +2*m,
          *pa12 = A +3*m, *pa13 = A +4*m, *pa23 = A +5*m;
    float *pap1 = Ap,     *pap2 = Ap+m,   *pap3 = Ap+2*m;
    float *pp1  = p,      *pp2  = p +m,   *pp3  = p +2*m;

    for(i=0; i<m; i++)
    {
        pap1[i] += pa11[i]*pp1[i] + pa12[i]*pp2[i] + pa13[i]*pp3[i];
        pap2[i] += pa12[i]*pp1[i] + pa22[i]*pp2[i] + pa23[i]*pp3[i];
        pap3[i] += pa13[i]*pp1[i] + pa23[i]*pp2[i] + pa33[i]*pp3[i];
    }
}

static double sumsq_le(int dm[], float a[], float b[], double s[], float u[])
{
    double ss = 0.0;
    int k;
    double mu = s[3], lam = s[4], id = s[5];
    double wx0, wx1, wx2, wy0, wy1, wy2, wz0, wz1, wz2, wxy, wxz, wyz;

    wx0 = 2*mu*(s[1]*s[1]+s[2]*s[2])+(4*mu+2*lam)*s[0]*s[0] + id;
    wy0 = 2*mu*(s[0]*s[0]+s[2]*s[2])+(4*mu+2*lam)*s[1]*s[1] + id;
    wz0 = 2*mu*(s[0]*s[0]+s[1]*s[1])+(4*mu+2*lam)*s[2]*s[2] + id;

    wx1 = -(2*mu+lam)*s[0]*s[0];
    wy1 = -(2*mu+lam)*s[1]*s[1];
    wz1 = -(2*mu+lam)*s[2]*s[2];

    wx2 = -mu*s[0]*s[0];
    wy2 = -mu*s[1]*s[1];
    wz2 = -mu*s[2]*s[2];

    wxy = 0.25*(lam+mu)*s[0]*s[1];
    wxz = 0.25*(lam+mu)*s[0]*s[2];
    wyz = 0.25*(lam+mu)*s[1]*s[2];

    for(k=0; k<dm[2]; k++)
    {
        int j, km1,kp1;
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *paxy, *payy, *paxz, *payz, *pazz;
            int i, jm1,jp1;

            pux  = u+dm[0]*(j+dm[1]*k);
            puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
            puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
            pbx  = b+dm[0]*(j+dm[1]*k);
            pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
            pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));
            paxx = a+dm[0]*(j+dm[1]*k);
            payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
            pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
            paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
            paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
            payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                int im1,ip1;
                float *px = pux+i, *py = puy+i, *pz = puz+i;
                double tmp;

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                tmp = (wx0+paxx[i])*px[0] + paxy[i]*py[0] + paxz[i]*pz[0]
                     + wx1*(px[ip1] + px[im1])
                     + wy2*(px[jp1] + px[jm1])
                     + wz2*(px[kp1] + px[km1])
                     + wxy*(py[ip1+jm1] - py[ip1+jp1] - py[im1+jm1] + py[im1+jp1])
                     + wxz*(pz[ip1+km1] - pz[ip1+kp1] - pz[im1+km1] + pz[im1+kp1])
                     - pbx[i];
                ss += tmp*tmp;

                tmp = (wy0+payy[i])*py[0] + paxy[i]*px[0] + payz[i]*pz[0]
                     + wx2*(py[ip1] + py[im1])
                     + wy1*(py[jp1] + py[jm1])
                     + wz2*(py[kp1] + py[km1])
                     + wxy*(px[jp1+im1] - px[jp1+ip1] - px[jm1+im1] + px[jm1+ip1])
                     + wyz*(pz[jp1+km1] - pz[jp1+kp1] - pz[jm1+km1] + pz[jm1+kp1])
                     - pby[i];
                ss += tmp*tmp;

                tmp = (wz0+pazz[i])*pz[0] + paxz[i]*px[0] + payz[i]*py[0]
                     + wx2*(pz[ip1] + pz[im1])
                     + wy2*(pz[jp1] + pz[jm1])
                     + wz1*(pz[kp1] + pz[km1])
                     + wxz*(px[kp1+im1] - px[kp1+ip1] - px[km1+im1] + px[km1+ip1])
                     + wyz*(py[kp1+jm1] - py[kp1+jp1] - py[km1+jm1] + py[km1+jp1])
                     - pbz[i];
                ss += tmp*tmp;

            }
        }
    }
    return(ss);
}

void LtLf_le(int dm[], float f[], double s[], float g[])
{
    int k;
    double mu = s[3], lam = s[4], id = s[5];
    double wx0, wx1, wx2, wy0, wy1, wy2, wz0, wz1, wz2, wxy, wxz, wyz;

    wx0 = 2*mu*(s[1]*s[1]+s[2]*s[2])+(4*mu+2*lam)*s[0]*s[0] + id;
    wy0 = 2*mu*(s[0]*s[0]+s[2]*s[2])+(4*mu+2*lam)*s[1]*s[1] + id;
    wz0 = 2*mu*(s[0]*s[0]+s[1]*s[1])+(4*mu+2*lam)*s[2]*s[2] + id;

    wx1 = -(2*mu+lam)*s[0]*s[0];
    wy1 = -(2*mu+lam)*s[1]*s[1];
    wz1 = -(2*mu+lam)*s[2]*s[2];

    wx2 = -mu*s[0]*s[0];
    wy2 = -mu*s[1]*s[1];
    wz2 = -mu*s[2]*s[2];

    wxy = 0.25*(lam+mu)*s[0]*s[1];
    wxz = 0.25*(lam+mu)*s[0]*s[2];
    wyz = 0.25*(lam+mu)*s[1]*s[2];

    for(k=0; k<dm[2]; k++)
    {
        int j, km1,kp1;
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            int i, jm1,jp1;
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
                int im1,ip1;
                float *px = &pfx[i], *py = &pfy[i], *pz = &pfz[i];

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                pgx[i] = wx0*px[0]
                       + wx1*(px[ip1] + px[im1])
                       + wy2*(px[jp1] + px[jm1])
                       + wz2*(px[kp1] + px[km1])
                       + wxy*(py[ip1+jm1] - py[ip1+jp1] - py[im1+jm1] + py[im1+jp1])
                       + wxz*(pz[ip1+km1] - pz[ip1+kp1] - pz[im1+km1] + pz[im1+kp1]);

                pgy[i] = wy0*py[0]
                       + wx2*(py[ip1] + py[im1])
                       + wy1*(py[jp1] + py[jm1])
                       + wz2*(py[kp1] + py[km1])
                       + wxy*(px[jp1+im1] - px[jp1+ip1] - px[jm1+im1] + px[jm1+ip1])
                       + wyz*(pz[jp1+km1] - pz[jp1+kp1] - pz[jm1+km1] + pz[jm1+kp1]);

                pgz[i] = wz0*pz[0]
                       + wx2*(pz[ip1] + pz[im1])
                       + wy2*(pz[jp1] + pz[jm1])
                       + wz1*(pz[kp1] + pz[km1])
                       + wxz*(px[kp1+im1] - px[kp1+ip1] - px[km1+im1] + px[km1+ip1])
                       + wyz*(py[kp1+jm1] - py[kp1+jp1] - py[km1+jm1] + py[km1+jp1]);
            }
        }
    }
}

static void relax_le(int dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double regx, regy, regz;
    double mu = s[3], lam = s[4], id = s[5];
    double wx0, wx1, wx2, wy0, wy1, wy2, wz0, wz1, wz2, wxy, wxz, wyz;

    wx0 = 2*mu*(s[1]*s[1]+s[2]*s[2])+(4*mu+2*lam)*s[0]*s[0] + id;
    wy0 = 2*mu*(s[0]*s[0]+s[2]*s[2])+(4*mu+2*lam)*s[1]*s[1] + id;
    wz0 = 2*mu*(s[0]*s[0]+s[1]*s[1])+(4*mu+2*lam)*s[2]*s[2] + id;

    wx1 = -(2*mu+lam)*s[0]*s[0];
    wy1 = -(2*mu+lam)*s[1]*s[1];
    wz1 = -(2*mu+lam)*s[2]*s[2];

    wx2 = -mu*s[0]*s[0];
    wy2 = -mu*s[1]*s[1];
    wz2 = -mu*s[2]*s[2];

    wxy = 0.25*(lam+mu)*s[0]*s[1];
    wxz = 0.25*(lam+mu)*s[0]*s[2];
    wyz = 0.25*(lam+mu)*s[1]*s[2];

    /* For stability in Gauss-Seidel relaxation, the magnitude of the diagonal element must
       exceed the sum of the magnitudes of the off diagonal elements of each column or row
       (see e.g. http://www.mathpages.com/home/kmath175/kmath175.htm).
       This should stabilise the relaxation, providing the second derives are positive definite. */
    regx = (4.0*(wxy+wxz)-2.0*(wx1+wy2+wz2)) - wx0; if (regx<0.0) regx = 0.0;
    regy = (4.0*(wxy+wyz)-2.0*(wx2+wy1+wz2)) - wy0; if (regy<0.0) regy = 0.0;
    regz = (4.0*(wxz+wyz)-2.0*(wx2+wy2+wz1)) - wz0; if (regz<0.0) regz = 0.0;

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("%dx%dx%d (%g,%g,%g): ", dm[0],dm[1],dm[2], regx,regy,regz);
#   endif

    for(it=0; it<8*nit; it++)
    {
        int k;
        /* double ss = 0.0; */
        for(k=it&1; k<dm[2]; k+=2)
        {
            int j, km1, kp1;
            km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

            for(j=(it>>1)&1; j<dm[1]; j+=2)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
                int i, jm1,jp1;

                pux  = u+dm[0]*(j+dm[1]* k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]* k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxx = a+dm[0]*(j+dm[1]* k);
                payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));

                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

                for(i=(it>>2)&1; i<dm[0]; i+=2)
                {
                    int im1,ip1;
                    double sux, suy, suz, axx, ayy, azz, axy, axz, ayz, idt;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;

                    sux = pbx[i] - ((wx0+paxx[i])*px[0] + paxy[i]*py[0] + paxz[i]*pz[0]
                                   + wx1*(px[ip1] + px[im1])
                                   + wy2*(px[jp1] + px[jm1])
                                   + wz2*(px[kp1] + px[km1])
                                   + wxy*(py[ip1+jm1] - py[ip1+jp1] - py[im1+jm1] + py[im1+jp1])
                                   + wxz*(pz[ip1+km1] - pz[ip1+kp1] - pz[im1+km1] + pz[im1+kp1]));

                    suy = pby[i] - ((wy0+payy[i])*py[0] + paxy[i]*px[0] + payz[i]*pz[0]
                                   + wx2*(py[ip1] + py[im1])
                                   + wy1*(py[jp1] + py[jm1])
                                   + wz2*(py[kp1] + py[km1])
                                   + wxy*(px[jp1+im1] - px[jp1+ip1] - px[jm1+im1] + px[jm1+ip1])
                                   + wyz*(pz[jp1+km1] - pz[jp1+kp1] - pz[jm1+km1] + pz[jm1+kp1]));

                    suz = pbz[i] - ((wz0+pazz[i])*pz[0] + paxz[i]*px[0] + payz[i]*py[0]
                                   + wx2*(pz[ip1] + pz[im1])
                                   + wy2*(pz[jp1] + pz[jm1])
                                   + wz1*(pz[kp1] + pz[km1])
                                   + wxz*(px[kp1+im1] - px[kp1+ip1] - px[km1+im1] + px[km1+ip1])
                                   + wyz*(py[kp1+jm1] - py[kp1+jp1] - py[km1+jm1] + py[km1+jp1]));

                    /* ss  += sux*sux + suy*suy + suz*suz; */

                    axx  = paxx[i] + wx0+regx;
                    ayy  = payy[i] + wy0+regy;
                    azz  = pazz[i] + wz0+regz;
                    axy  = paxy[i];
                    axz  = paxz[i];
                    ayz  = payz[i];
                    idt  = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);

                    *px += idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                    *py += idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                    *pz += idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));

                }
            }
        }

#       ifdef VERBOSE
            printf(" %g", sumsq_le(dm, a, b, s, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}


static void Atimesp_le(int dm[], float A[], double param[], float p[], float Ap[])
{
    LtLf_le(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}


static double sumsq_me(int dm[], float a[], float b[], double s[], float u[])
{
    double w000, w001, w010, w100;
    double ss = 0.0;
    int i, j, k;

    w000 = s[3]*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + s[5];
    w001 = s[3]*(-s[2]*s[2]);
    w010 = s[3]*(-s[1]*s[1]);
    w100 = s[3]*(-s[0]*s[0]);

    for(k=0; k<dm[2]; k++)
    {
        int km1,kp1;
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *paxy, *payy, *paxz, *payz, *pazz;
            int jm1,jp1,im1,ip1;

            pux  = u+dm[0]*(j+dm[1]*k);
            puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
            puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
            pbx  = b+dm[0]*(j+dm[1]*k);
            pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
            pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));
            paxx = a+dm[0]*(j+dm[1]*k);
            payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
            pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
            paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
            paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
            payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                float *px = &pux[i], *py = &puy[i], *pz = &puz[i];
                double tmp;

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                tmp = (w000+paxx[i])*px[0] + paxy[i]*py[0] + paxz[i]*pz[0]
                     + w001*(px[km1] + px[kp1])
                     + w010*(px[jm1] + px[jp1])
                     + w100*(px[im1] + px[ip1])
                     - pbx[i];
                ss += tmp*tmp;

                tmp = (w000+payy[i])*py[0] + paxy[i]*px[0] + payz[i]*pz[0]
                     + w001*(py[km1] + py[kp1])
                     + w010*(py[jm1] + py[jp1])
                     + w100*(py[im1] + py[ip1])
                     - pby[i];
                ss += tmp*tmp;

                tmp = (w000+pazz[i])*pz[0] + paxz[i]*px[0] + payz[i]*py[0]
                     + w001*(pz[km1] + pz[kp1])
                     + w010*(pz[jm1] + pz[jp1])
                     + w100*(pz[im1] + pz[ip1])
                     - pbz[i];
                ss += tmp*tmp;
            }
        }
    }
    return(ss);
}

void LtLf_me(int dm[], float f[], double s[], float g[])
{
    int i, j, k, km1,kp1, jm1,jp1, im1,ip1;
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

                pgx[i] = w000*px[0] + w001*(px[km1] + px[kp1]) + w010*(px[jm1] + px[jp1]) + w100*(px[im1] + px[ip1]);
                pgy[i] = w000*py[0] + w001*(py[km1] + py[kp1]) + w010*(py[jm1] + py[jp1]) + w100*(py[im1] + py[ip1]);
                pgz[i] = w000*pz[0] + w001*(pz[km1] + pz[kp1]) + w010*(pz[jm1] + pz[jp1]) + w100*(pz[im1] + pz[ip1]);
            }
        }
    }
}

static void relax_me(int dm[], float a[], float b[], double s[], int nit, float u[])
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
        int k, kstart;
        int j, jstart;
        int i, istart;

#       ifdef VERBOSE
            printf(" %g", sumsq_me(dm, a, b, s, u));
#       endif

        kstart = it%2;
        for(k=0; k<dm[2]; k++)
        {
            int km1, kp1;
            km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

            jstart = (kstart == (k%2));
            for(j=0; j<dm[1]; j++)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *paxy, *payy, *paxz, *payz, *pazz;
                int jm1,jp1, im1,ip1;

                pux  = u+dm[0]*(j+dm[1]*k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]*k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxx = a+dm[0]*(j+dm[1]*k);
                payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));

                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

                istart = (jstart == (j%2));

                for(i=istart; i<dm[0]; i+=2)
                {
                    double sux, suy, suz, axx, ayy, azz, axy, axz, ayz, idt;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;

                    sux = pbx[i]-(w001*(px[km1] + px[kp1]) + w010*(px[jm1] + px[jp1]) + w100*(px[im1] + px[ip1]));
                    suy = pby[i]-(w001*(py[km1] + py[kp1]) + w010*(py[jm1] + py[jp1]) + w100*(py[im1] + py[ip1]));
                    suz = pbz[i]-(w001*(pz[km1] + pz[kp1]) + w010*(pz[jm1] + pz[jp1]) + w100*(pz[im1] + pz[ip1]));
                    /*
                       syms axx ayy azz axy axz ayz sux suy suz
                       A = [axx axy axz; axy ayy ayz; axz ayz azz];
                       su = [sux ; suy; suz]
                       simplify(inv(A)*su)
                    */
                    axx = paxx[i] + w000;
                    ayy = payy[i] + w000;
                    azz = pazz[i] + w000;
                    axy = paxy[i];
                    axz = paxz[i];
                    ayz = payz[i];
                    idt = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);
                    *px = idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                    *py = idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                    *pz = idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                }
            }
        }
    }
#   ifdef VERBOSE
        printf(" %g\n", sumsq_me(dm, a, b, s, u));
#   endif
}

static void Atimesp_me(int dm[], float A[], double param[], float p[], float Ap[])
{
    LtLf_me(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

static double sumsq_be(int dm[], float a[], float b[], double s[], float u[])
{
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double ss = 0.0;
    int k;

    w000 = s[3]*(6*(s[0]*s[0]*s[0]*s[0]+s[1]*s[1]*s[1]*s[1]+s[2]*s[2]*s[2]*s[2])
                +8*(s[0]*s[0]*s[1]*s[1]+s[0]*s[0]*s[2]*s[2]+s[1]*s[1]*s[2]*s[2])
                +4*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])*s[4]  + s[4]*s[4]) + s[5];
    w100 = s[3]*(-2*s[0]*s[0]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w200 = s[3]*s[0]*s[0]*s[0]*s[0];
    w010 = s[3]*(-2*s[1]*s[1]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w110 = s[3]*2*s[0]*s[0]*s[1]*s[1];
    w020 = s[3]*s[1]*s[1]*s[1]*s[1];
    w001 = s[3]*(-2*s[2]*s[2]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w101 = s[3]*2*s[0]*s[0]*s[2]*s[2];
    w011 = s[3]*2*s[1]*s[1]*s[2]*s[2];
    w002 = s[3]*s[2]*s[2]*s[2]*s[2];

    for(k=0; k<dm[2]; k++)
    {
        int j, km2,km1,kp1,kp2;
        km2 = (BOUND(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (BOUND(k+2,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
            int i, jm2,jm1,jp1,jp2;

            pux  = u+dm[0]*(j+dm[1]*k);
            puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
            puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
            pbx  = b+dm[0]*(j+dm[1]*k);
            pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
            pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));
            paxx = a+dm[0]*(j+dm[1]*k);
            payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
            pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
            paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
            paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
            payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));

            jm2 = (BOUND(j-2,dm[1])-j)*dm[0];
            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];
            jp2 = (BOUND(j+2,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                int im2,im1,ip1,ip2;
                float *px = pux+i, *py = puy+i, *pz = puz+i;
                double tmp;

                im2 = BOUND(i-2,dm[0])-i;
                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;
                ip2 = BOUND(i+2,dm[0])-i;

                tmp = (w000+paxx[i])*px[0] + paxy[i]*py[0] + paxz[i]*pz[0]
                     + w010*(px[    jm1    ] + px[    jp1    ])
                     + w020*(px[    jm2    ] + px[    jp2    ])
                     + w100*(px[im1        ] + px[ip1        ])
                     + w110*(px[im1+jm1    ] + px[ip1+jm1    ] + px[im1+jp1    ] + px[ip1+jp1    ])
                     + w200*(px[im2        ] + px[ip2        ])
                     + w001*(px[        km1] + px[        kp1])
                     + w101*(px[im1    +km1] + px[ip1    +km1] + px[im1    +kp1] + px[ip1    +kp1])
                     + w011*(px[    jm1+km1] + px[    jp1+km1] + px[    jm1+kp1] + px[    jp1+kp1])
                     + w002*(px[        km2] + px[        kp2])
                     - pbx[i];
                ss += tmp*tmp;

                tmp = paxy[i]*px[0] + (w000+payy[i])*py[0] + payz[i]*pz[0]
                     + w010*(py[    jm1    ] + py[    jp1    ])
                     + w020*(py[    jm2    ] + py[    jp2    ])
                     + w100*(py[im1        ] + py[ip1        ])
                     + w110*(py[im1+jm1    ] + py[ip1+jm1    ] + py[im1+jp1    ] + py[ip1+jp1    ])
                     + w200*(py[im2        ] + py[ip2        ])
                     + w001*(py[        km1] + py[        kp1])
                     + w101*(py[im1    +km1] + py[ip1    +km1] + py[im1    +kp1] + py[ip1    +kp1])
                     + w011*(py[    jm1+km1] + py[    jp1+km1] + py[    jm1+kp1] + py[    jp1+kp1])
                     + w002*(py[        km2] + py[        kp2])
                     - pby[i];
                ss += tmp*tmp;

                tmp = paxz[i]*px[0] + payz[i]*py[0] + (w000+pazz[i])*pz[0]
                     + w010*(pz[    jm1    ] + pz[    jp1    ])
                     + w020*(pz[    jm2    ] + pz[    jp2    ])
                     + w100*(pz[im1        ] + pz[ip1        ])
                     + w110*(pz[im1+jm1    ] + pz[ip1+jm1    ] + pz[im1+jp1    ] + pz[ip1+jp1    ])
                     + w200*(pz[im2        ] + pz[ip2        ])
                     + w001*(pz[        km1] + pz[        kp1])
                     + w101*(pz[im1    +km1] + pz[ip1    +km1] + pz[im1    +kp1] + pz[ip1    +kp1])
                     + w011*(pz[    jm1+km1] + pz[    jp1+km1] + pz[    jm1+kp1] + pz[    jp1+kp1])
                     + w002*(pz[        km2] + pz[        kp2])
                     - pbz[i];
                ss += tmp*tmp;
            }
        }
    }
    return(ss);
}

void LtLf_be(int dm[], float f[], double s[], float g[])
{
    int k;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    /*
        syms s1 s2 s3
        zz = sym(zeros(3,3));
        K1 = cat(3,zz,[0 -s1*s1 0; 0 2*s1*s1 0; 0 -s1*s1 0],zz);
        K2 = cat(3,zz,[0 0 0; -s2*s2 2*s2*s2 -s2*s2; 0 0 0],zz);
        K3 = sym(zeros(3,3,3));
        K3(2,2,1) = -s3*s3;
        K3(2,2,2) = 2*s3*s3;
        K3(2,2,3) = -s3*s3;

        K  = K1+K2+K3;
        % L  = convn(K,K)
        L  = sym(zeros(5,5,5));
        for i=1:3,
            for j=1:3,
                for k=1:3,
                    L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) = L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) + K(i,j,k)*K;
                end;
            end;
        end;
        disp(L(3:end,3:end,3:end))
    */

    w000 = s[3]*(6*(s[0]*s[0]*s[0]*s[0]+s[1]*s[1]*s[1]*s[1]+s[2]*s[2]*s[2]*s[2])
                +8*(s[0]*s[0]*s[1]*s[1]+s[0]*s[0]*s[2]*s[2]+s[1]*s[1]*s[2]*s[2])
                +4*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])*s[4]  + s[4]*s[4]) + s[5];
    w100 = s[3]*(-2*s[0]*s[0]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w200 = s[3]*s[0]*s[0]*s[0]*s[0];
    w010 = s[3]*(-2*s[1]*s[1]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w110 = s[3]*2*s[0]*s[0]*s[1]*s[1];
    w020 = s[3]*s[1]*s[1]*s[1]*s[1];
    w001 = s[3]*(-2*s[2]*s[2]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w101 = s[3]*2*s[0]*s[0]*s[2]*s[2];
    w011 = s[3]*2*s[1]*s[1]*s[2]*s[2];
    w002 = s[3]*s[2]*s[2]*s[2]*s[2];

    for(k=0; k<dm[2]; k++)
    {
        int j, km2,km1,kp1,kp2;
        km2 = (BOUND(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (BOUND(k+2,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            int i, jm2,jm1,jp1,jp2;
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
                int im2,im1,ip1,ip2;
                float *px = &pfx[i], *py = &pfy[i], *pz = &pfz[i];

                im2 = BOUND(i-2,dm[0])-i;
                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;
                ip2 = BOUND(i+2,dm[0])-i;

                pgx[i] = w000* px[0]
                       + w010*(px[    jm1    ] + px[    jp1    ])
                       + w020*(px[    jm2    ] + px[    jp2    ])
                       + w100*(px[im1        ] + px[ip1        ])
                       + w110*(px[im1+jm1    ] + px[ip1+jm1    ] + px[im1+jp1    ] + px[ip1+jp1    ])
                       + w200*(px[im2        ] + px[ip2        ])
                       + w001*(px[        km1] + px[        kp1])
                       + w101*(px[im1    +km1] + px[ip1    +km1] + px[im1    +kp1] + px[ip1    +kp1])
                       + w011*(px[    jm1+km1] + px[    jp1+km1] + px[    jm1+kp1] + px[    jp1+kp1])
                       + w002*(px[        km2] + px[        kp2]);
                pgy[i] = w000* py[0]
                       + w010*(py[    jm1    ] + py[    jp1    ])
                       + w020*(py[    jm2    ] + py[    jp2    ])
                       + w100*(py[im1        ] + py[ip1        ])
                       + w110*(py[im1+jm1    ] + py[ip1+jm1    ] + py[im1+jp1    ] + py[ip1+jp1    ])
                       + w200*(py[im2        ] + py[ip2        ])
                       + w001*(py[        km1] + py[        kp1])
                       + w101*(py[im1    +km1] + py[ip1    +km1] + py[im1    +kp1] + py[ip1    +kp1])
                       + w011*(py[    jm1+km1] + py[    jp1+km1] + py[    jm1+kp1] + py[    jp1+kp1])
                       + w002*(py[        km2] + py[        kp2]);
                pgz[i] = w000* pz[0]
                       + w010*(pz[    jm1    ] + pz[    jp1    ])
                       + w020*(pz[    jm2    ] + pz[    jp2    ])
                       + w100*(pz[im1        ] + pz[ip1        ])
                       + w110*(pz[im1+jm1    ] + pz[ip1+jm1    ] + pz[im1+jp1    ] + pz[ip1+jp1    ])
                       + w200*(pz[im2        ] + pz[ip2        ])
                       + w001*(pz[        km1] + pz[        kp1])
                       + w101*(pz[im1    +km1] + pz[ip1    +km1] + pz[im1    +kp1] + pz[ip1    +kp1])
                       + w011*(pz[    jm1+km1] + pz[    jp1+km1] + pz[    jm1+kp1] + pz[    jp1+kp1])
                       + w002*(pz[        km2] + pz[        kp2]);
            }
        }
    }
}

static void relax_be(int dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double reg;

    w000 = s[3]*(6*(s[0]*s[0]*s[0]*s[0]+s[1]*s[1]*s[1]*s[1]+s[2]*s[2]*s[2]*s[2])
                +8*(s[0]*s[0]*s[1]*s[1]+s[0]*s[0]*s[2]*s[2]+s[1]*s[1]*s[2]*s[2])
                +4*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])*s[4]  + s[4]*s[4]) + s[5];
    w100 = s[3]*(-2*s[0]*s[0]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w200 = s[3]*s[0]*s[0]*s[0]*s[0];
    w010 = s[3]*(-2*s[1]*s[1]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w110 = s[3]*2*s[0]*s[0]*s[1]*s[1];
    w020 = s[3]*s[1]*s[1]*s[1]*s[1];
    w001 = s[3]*(-2*s[2]*s[2]*(2*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])+s[4]));
    w101 = s[3]*2*s[0]*s[0]*s[2]*s[2];
    w011 = s[3]*2*s[1]*s[1]*s[2]*s[2];
    w002 = s[3]*s[2]*s[2]*s[2]*s[2];

    /* For stability in Gauss-Seidel relaxation, the magnitude of the diagonal element must
       exceed the sum of the magnitudes of the off diagonal elements of each column or row
       (see e.g. http://www.mathpages.com/home/kmath175/kmath175.htm).
       This is an attempt to stabilise the relaxation
       NOTE: It still isn't completely stable and can diverge as well as converge. */
    reg = (2.0*(w200+w020+w002)-2.0*(w100+w010+w001)+4.0*(w110+w011+w101)) - w000;
    if (reg<0.0) reg = 0.0;

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("%dx%dx%d (%g): ", dm[0],dm[1],dm[2],reg);
#   endif

    for(it=0; it<8*nit; it++)
    {
        int k, kstart,kend,kskip;
        int j, jstart,jend,jskip;
        int i, istart,iend,iskip;
        /* double ss = 0.0; */

        if ((it/2/2)%2)
        {
            kstart = 0;
            kend   = dm[2];
            kskip  = 1;
        }
        else
        {
            kstart = dm[2]-1;
            kend   = -1;
            kskip  = -1;
        }
        if ((it/2)%2)
        {
            jstart = 0;
            jend   = dm[1];
            jskip  = 1;
        }
        else
        {
            jstart = dm[1]-1;
            jend   = -1;
            jskip  = -1;
        }
        if (it%2)
        {
            istart = 0;
            iend   = dm[0];
            iskip  = 1;
        }
        else
        {
            istart = dm[0]-1;
            iend   = -1;
            iskip  = -1;
        }

        for(k=kstart; k!=kend; k+=kskip)
        {
            int km2, km1, kp1, kp2;
            km2 = (BOUND(k-2,dm[2])-k)*dm[0]*dm[1];
            km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];
            kp2 = (BOUND(k+2,dm[2])-k)*dm[0]*dm[1];

            for(j=jstart; j!=jend; j+=jskip)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
                int jm2,jm1,jp1,jp2;

                pux  = u+dm[0]*(j+dm[1]* k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]* k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxx = a+dm[0]*(j+dm[1]* k);
                payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));

                jm2 = (BOUND(j-2,dm[1])-j)*dm[0];
                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];
                jp2 = (BOUND(j+2,dm[1])-j)*dm[0];

                for(i=istart; i!=iend; i+=iskip)
                {
                    int im2,im1,ip1,ip2;
                    double sux, suy, suz, axx, ayy, azz, axy, axz, ayz, idt;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im2 = BOUND(i-2,dm[0])-i;
                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;
                    ip2 = BOUND(i+2,dm[0])-i;

                    sux = pbx[i] - ((w000+paxx[i])*px[0] + paxy[i]*py[0] + paxz[i]*pz[0]
                                  + w010*(px[    jm1    ] + px[    jp1    ])
                                  + w020*(px[    jm2    ] + px[    jp2    ])
                                  + w100*(px[im1        ] + px[ip1        ])
                                  + w110*(px[im1+jm1    ] + px[ip1+jm1    ] + px[im1+jp1    ] + px[ip1+jp1    ])
                                  + w200*(px[im2        ] + px[ip2        ])
                                  + w001*(px[        km1] + px[        kp1])
                                  + w101*(px[im1    +km1] + px[ip1    +km1] + px[im1    +kp1] + px[ip1    +kp1])
                                  + w011*(px[    jm1+km1] + px[    jp1+km1] + px[    jm1+kp1] + px[    jp1+kp1])
                                  + w002*(px[        km2] + px[        kp2]));

                    suy = pby[i] - (paxy[i]*px[0] + (w000+payy[i])*py[0] + payz[i]*pz[0]
                                  + w010*(py[    jm1    ] + py[    jp1    ])
                                  + w020*(py[    jm2    ] + py[    jp2    ])
                                  + w100*(py[im1        ] + py[ip1        ])
                                  + w110*(py[im1+jm1    ] + py[ip1+jm1    ] + py[im1+jp1    ] + py[ip1+jp1    ])
                                  + w200*(py[im2        ] + py[ip2        ])
                                  + w001*(py[        km1] + py[        kp1])
                                  + w101*(py[im1    +km1] + py[ip1    +km1] + py[im1    +kp1] + py[ip1    +kp1])
                                  + w011*(py[    jm1+km1] + py[    jp1+km1] + py[    jm1+kp1] + py[    jp1+kp1])
                                  + w002*(py[        km2] + py[        kp2]));

                    suz = pbz[i] - (paxz[i]*px[0] + payz[i]*py[0] + (w000+pazz[i])*pz[0]
                                  + w010*(pz[    jm1    ] + pz[    jp1    ])
                                  + w020*(pz[    jm2    ] + pz[    jp2    ])
                                  + w100*(pz[im1        ] + pz[ip1        ])
                                  + w110*(pz[im1+jm1    ] + pz[ip1+jm1    ] + pz[im1+jp1    ] + pz[ip1+jp1    ])
                                  + w200*(pz[im2        ] + pz[ip2        ])
                                  + w001*(pz[        km1] + pz[        kp1])
                                  + w101*(pz[im1    +km1] + pz[ip1    +km1] + pz[im1    +kp1] + pz[ip1    +kp1])
                                  + w011*(pz[    jm1+km1] + pz[    jp1+km1] + pz[    jm1+kp1] + pz[    jp1+kp1])
                                  + w002*(pz[        km2] + pz[        kp2]));
                    /* ss  += sux*sux + suy*suy + suz*suz; */

                    axx  = paxx[i] + w000+reg;
                    ayy  = payy[i] + w000+reg;
                    azz  = pazz[i] + w000+reg;
                    axy  = paxy[i];
                    axz  = paxz[i];
                    ayz  = payz[i];
                    idt  = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);
                    *px += idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                    *py += idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                    *pz += idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                }
            }
        }
#       ifdef VERBOSE
            printf(" %g", sumsq_be(dm, a, b, s, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}


static void Atimesp_be(int dm[], float A[], double param[], float p[], float Ap[])
{
    LtLf_be(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

/*
syms a0 a1 a2 a3 a4 a5 b0 b1 b2
A = [a0 a3 a4; a3 a1 a5; a4 a5 a2];
b = [b0 b1 b2].';
A\b
*/
static void solve33(float a[], float b[], double t, float u[])
{
    double dt;
    double a0 = a[0]+t, a1 = a[1]+t, a2 = a[2]+t;

    dt  = a0*a2*a1-a0*a[5]*a[5]-a1*a[4]*a[4]-a2*a[3]*a[3]+2*a[5]*a[3]*a[4];
    if (dt<1e-64*a0*a1*a2)
        dt = 1e-64*a0*a1*a2;
    u[0] = (b[0]*(a2*a1-a[5]*a[5])+b[1]*(a[4]*a[5]-a[3]*a2)+b[2]*(a[3]*a[5]-a[4]*a1))/dt;
    u[1] = (b[0]*(a[5]*a[4]-a[3]*a2)+b[1]*(a0*a2-a[4]*a[4])+b[2]*(a[4]*a[3]-a0*a[5]))/dt;
    u[2] = (b[0]*(a[5]*a[3]-a[4]*a1)+b[1]*(a[4]*a[3]-a0*a[5])+b[2]*(a0*a1-a[3]*a[3]))/dt;
}

static double dotprod(int m, float a[], float b[])
{
    int i;
    double dp = 0.0;
    for(i=0; i<m; i++)
        dp += a[i]*b[i];
    return(dp);
}

static void addscaled(int m, float a[], float b[], double s)
{
    int i;
    for(i=0; i<m; i++)
        a[i] += s*b[i];
}

float norm(int m, float a[])
{
    int i;
    double dp = 0.0;
    for(i=0; i<m; i++)
        dp += a[i]*a[i];
    return(sqrt(dp));
}


/*
% Solve A*x = b by the conjugate gradient method
% See Gilbert, Moler & Schreiber (1991)
% Sparse Matrices in Matlab: Design and Implementation
% SIAM Journal on Matrix Analysis and Applications
% http://citeseer.ist.psu.edu/gilbert91sparse.html

if nargin<3, tol = 1e-4; end;
if nargin<4, nit = 1000; end;

x    = zeros(size(b));
r    = b;
rtr  = r'*r;
p    = zeros(size(b));
beta = 0;
it   = 0;
while norm(r) > tol*norm(b),
    p      = r + beta*p;
    Ap     = A*p;
    alpha  = rtr/(p'*Ap);
    x      = x + alpha*p;
    r      = r - alpha*Ap;
    rtrold = rtr;
    rtr    = r'*r;
    beta   = rtr/rtrold;

    it = it+1;
    if it>nit, break; end;
end;
*/

void cgs3(int dm[], float A[], float b[], int rtype, double param[], double tol, int nit,
             float x[], float r[], float p[], float Ap[])
{
    int i, m = dm[0]*dm[1]*dm[2]*3, it;
    double rtr, nb, rtrold, alpha, beta;
    void (*Atimesp)();

    /* printf("\n **** %dx%d ****\n",dm[0],dm[1]); */
    if (rtype == 0)
        Atimesp = Atimesp_le;
    else if (rtype == 1)
        Atimesp = Atimesp_me;
    else
        Atimesp = Atimesp_be;

    nb      = tol*norm(m,b);

#   ifdef NEVER
        /* Assuming starting estimates of zeros */
        /* x    = zeros(size(b)); */
        for(i=0; i<m;i++)
            x[i] = 0.0;

        /* r    = b; */
        for(i=0; i<m;i++)
            r[i] = b[i];
#   else
        /* Assume starting estimates are passed as arguments */
        /* r    = b-A*x; */
        Atimesp(dm, A, param, x, Ap);
        for(i=0; i<m;i++)
            r[i] = b[i]-Ap[i];
#   endif

    /* rtr  = r'*r; */
    rtr     = dotprod(m, r, r);

    /* p    = zeros(size(b)); */
    for(i=0; i<m;i++)
        p[i] = 0.0;

    /* beta = 0; */
    beta    = 0.0;

    /* for it=1:nit, */
    for(it=0; it<nit; it++)
    {
        /* if norm(r) < tol*norm(b), break; end; */
        if (norm(m,r) < nb)
            break;

        /* p      = r + beta*p; */
        for(i=0; i<m; i++)
            p[i]  = r[i] + beta*p[i];

        /* Ap     = A*p; */
        Atimesp(dm, A, param, p, Ap);

        /* alpha  = rtr/(p'*Ap); */
        alpha     = rtr/dotprod(m, p, Ap);

        /* x      = x + alpha*p; */
        addscaled(m, x, p, alpha);

        /* r      = r - alpha*Ap; */
        addscaled(m, r, Ap, -alpha);

        /* rtrold = rtr; */
        rtrold = rtr;

        /* rtr    = r'*r; */
        rtr       = dotprod(m, r, r);

        /* beta   = rtr/rtrold; */
        beta      = rtr/rtrold;

        /* printf("%d\t%g\t%g  %g %g\n",it, norm(m,r), nb/tol, alpha, beta); */
    /* end; */
    }
    /* printf("Done after %d iterations (%g, %g).\n",it, norm(m,r), nb); */
}

/*******************************************************/


static double wt2(double x)
{
        x = fabs(x);
        if (x < 0.5)
                return(0.75 - x*x);
        if (x < 1.5)
        {
                x = 1.5 - x;
                return(0.5*x*x);
        }
        return(0.0);
}

static void resized_plane(int na[], float *a,  int nc[], float *c, float *b)
{
    int i, j, o,om,op;
    double loc, s, w, wm, wp;
    float *ap, *bp, *cp;
    /* a - na[0]*na[1]
     * c - nc[0]*nc[1]
     * b - na[0]*nc[1]
     */

    s = (double)na[1]/(double)nc[1];
    for(j=0; j<nc[1]; j++)
    {
        loc = (j+0.5)*s-0.5;
        o   = floor(loc+0.5);
        om  = BOUND(o-1,na[1])*na[0];
        op  = BOUND(o+1,na[1])*na[0];
        w   = wt2( o   -loc);
        wp  = wt2((o+1)-loc);
        wm  = wt2((o-1)-loc);
        o  *= na[0];
        for(ap=a, bp=b+j*na[0], cp=ap+na[0]; ap<cp; ap++, bp++)
            *bp = wm*ap[om]+w*ap[o]+wp*ap[op];
    }
    s = (double)na[0]/(double)nc[0];
    for(i=0; i<nc[0]; i++)
    {
        loc = (i+0.5)*s-0.5;
        o   = floor(loc+0.5);
        om  = BOUND(o-1,na[0]);
        op  = BOUND(o+1,na[0]);
        w   = wt2( o   -loc);
        wp  = wt2((o+1)-loc);
        wm  = wt2((o-1)-loc);
        for(bp=b, cp=c+i, ap=bp+na[0]*nc[1]; bp<ap; bp+=na[0], cp+=nc[0])
            *cp = wm*bp[om]+w*bp[o]+wp*bp[op];
    }
}

void resize(int na[], float *a,  int nc[], float *c, float *b)
{
    int j, k, o=-999999,om,op, m, oo;
    double loc, s, w, wm, wp;
    float *bp, *cp, *pl[3];
    m     = nc[0]*nc[1];
    pl[0] = b;
    pl[1] = b + m;
    pl[2] = b + m*2;
    bp    = b + m*3;

    for(k=0; k<nc[2]; k++)
    {
        s      = (double)na[2]/(double)nc[2];
        loc    = (k+0.5)*s-0.5;
        oo     = o;
        o      = floor(loc+0.5);
        om     = BOUND(o-1,na[2]);
        op     = BOUND(o+1,na[2]);

        if (o==oo)
        {   /* do nothing */
        }
        else if (o==oo+1)
        {   /* Shift by 1 */
            float *tp;
            tp    = pl[0];
            pl[0] = pl[1];
            pl[1] = pl[2];
            pl[2] = tp;
            resized_plane(na, a+na[0]*na[1]*op,nc,pl[2],bp);
        }
        else if (o==oo+2)
        {   /* Shift by 2 */
            float *tp;
            tp    = pl[0];
            pl[0] = pl[2];
            pl[2] = tp;
            resized_plane(na, a+na[0]*na[1]*o ,nc,pl[1],bp);
            resized_plane(na, a+na[0]*na[1]*op,nc,pl[2],bp);
        }
        else
        {   /* Read everything */
            resized_plane(na, a+na[0]*na[1]*om,nc,pl[0],bp);
            resized_plane(na, a+na[0]*na[1]*o ,nc,pl[1],bp);
            resized_plane(na, a+na[0]*na[1]*op,nc,pl[2],bp);
        }
        w   = wt2( o   -loc);
        wp  = wt2((o+1)-loc);
        wm  = wt2((o-1)-loc);
        cp  = c+nc[0]*nc[1]*k;
        for(j=0; j<nc[0]*nc[1]; j++)
        {
            cp[j] = wm*pl[0][j]+w*pl[1][j]+wp*pl[2][j];
        }
    }
}

static void rescale(int n, float *a, double s)
{
    int i;
    for(i=0; i<n; i++)
        a[i] *= s;
}

static void restrict(int n,  int na[], float *a,  int nc[], float *c, float *b)
{
    int i;
    for(i=0; i<n; i++)
    {
        resize(na, a+i*na[0]*na[1]*na[2], nc, c+i*nc[0]*nc[1]*nc[2], b);
    }
}

static void prolong(int n,  int na[], float *a,  int nc[], float *c, float *b)
{
    int i;
    for(i=0; i<n; i++)
        resize(na, a+i*na[0]*na[1]*na[2], nc, c+i*nc[0]*nc[1]*nc[2], b);
}

static void zeros(int n, float *a)
{
    int i;
    for(i=0; i<n; i++)
        a[i] = 0.0;
}

static void copy(int n, float *a, float *b)
{
    int i;
    for(i=0; i<n; i++)
        b[i] = a[i];
}

static void addto(int n, float *a, float *b)
{
    int i;
    for(i=0; i<n; i++)
        a[i] += b[i];
}

int fmg3_scratchsize(int n0[])
{
    int    n[32][3], m[32], bs, j;
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
    return((3*n0[0]*n0[1]*n0[2] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 15*bs));
}

/*
    Full Multigrid solver.  See Numerical Recipes (second edition) for more
    information
*/
void fmg3(int n0[], float *a0, float *b0, int rtype, double param0[], int c, int nit,
          float *u0, float *scratch)
{
    int i, j, ng, bs;
     int n[32][3], m[32];
    float *bo[32], *a[32], *b[32], *u[32], *res, *rbuf;
    double param[32][6];
    void (*relax)(), (*Atimesp)();
    double (*sumsq)();

    if (rtype == 0)
    {
        relax   = relax_le;
        Atimesp = Atimesp_le;
        sumsq   = sumsq_le;
    }
    else if (rtype == 1)
    {
        relax   = relax_me;
        Atimesp = Atimesp_me;
        sumsq   = sumsq_me;
    }
    else
    {
        relax   = relax_be;
        Atimesp = Atimesp_be;
        sumsq   = sumsq_be;
    }

#   ifdef VERVOSE
        printf("start=%g\n", sumsq(n0, a0, b0, param[0], u0));
#   endif

    bo[0]   = b0;
    b[0]    = b0;
    u[0]    = u0;
    a[0]    = a0;
    n[0][0] = n0[0];
    n[0][1] = n0[1];
    n[0][2] = n0[2];
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
        m[j]    = n[j][0]*n[j][1]*n[j][2];
        ng ++;
        bs += m[j];
        if ((n[j][0]<2) && (n[j][1]<2) && (n[j][2]<2))
            break;
    }

    res    = scratch;
    rbuf   = scratch + 3*m[0];
    bo[1]  = scratch + 3*m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1];
    b[1]   = scratch + 3*m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 3*bs;
    u[1]   = scratch + 3*m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 6*bs;
    a[1]   = scratch + 3*m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 9*bs;

    for(j=2; j<ng; j++)
    {
        bo[j] = bo[j-1]+3*m[j-1];
        b[j]  =  b[j-1]+3*m[j-1];
        u[j]  =  u[j-1]+3*m[j-1];
        a[j]  =  a[j-1]+6*m[j-1];
    }

    for(j=1; j<ng; j++)
    {
        restrict(3,n[j-1],bo[j-1],n[j],bo[j],rbuf);
        restrict(6,n[j-1],a[j-1],n[j],a[j],rbuf);

        param[j][0] = param0[0]*(double)n[j][0]/n0[0];
        param[j][1] = param0[1]*(double)n[j][1]/n0[1];
        param[j][2] = param0[2]*(double)n[j][2]/n0[2];
        param[j][3] = param[0][3];
        param[j][4] = param[0][4];
        param[j][5] = param[0][5];
    }

    if (u[0][0]==0)
    {
     /* solve33(a[ng-1], bo[ng-1], param0[5], u[ng-1]); */
        relax(n[ng-1], a[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]);

        for(j=ng-2; j>=0; j--)
        {
            int jc;
            prolong(3,n[j+1],u[j+1],n[j],u[j],rbuf);
            if(j>0) copy(3*m[j],bo[j],b[j]);
            for(jc=0; jc<c; jc++)
            {
                int jj;
                for(jj=j; jj<ng-1; jj++)
                {
                    relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
                    Atimesp(n[jj], a[jj], param[jj], u[jj], res);
                    for(i=0; i<3*m[jj]; i++)
                        res[i] = b[jj][i] - res[i];
                    restrict(3,n[jj],res,n[jj+1],b[jj+1],rbuf);
                    zeros(3*m[jj+1],u[jj+1]);
                }
             /* solve33(a[ng-1], b[ng-1], param0[5], u[ng-1]); */
                relax(n[ng-1], a[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]);
                for(jj=ng-2; jj>=j; jj--)
                {
                    prolong(3,n[jj+1],u[jj+1],n[jj],res,rbuf);
                    addto(3*m[jj], u[jj], res);
                    relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
                }
            }
        }
    }
    else
    {
        int jc;
        for(j=1; j<ng; j++)
            restrict(3,n[j-1],u[j-1],n[j],u[j],rbuf);

        for(jc=0; jc<c; jc++)
        {
            int jj;
            for(jj=0; jj<ng-1; jj++)
            {
                relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
                Atimesp(n[jj], a[jj], param[jj], u[jj], res);
                for(i=0; i<3*m[jj]; i++)
                    res[i] = b[jj][i] - res[i];
                restrict(3,n[jj],res,n[jj+1],b[jj+1],rbuf);
                zeros(3*m[jj+1],u[jj+1]);
            }
         /* solve33(a[ng-1], b[ng-1], param0[5], u[ng-1]); */
            relax(n[ng-1], a[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]);

            for(jj=ng-2; jj>=0; jj--)
            {
                prolong(3,n[jj+1],u[jj+1],n[jj],res,rbuf);
                addto(3*m[jj], u[jj], res);
                relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
            }
        }
    }
    /* printf("end=%g\n", sumsq(n0, a0, b0, param[0], u0)); */
}

int fmg3_scratchsize_noa(int n0[])
{
    int    n[32][3], m[32], bs, j;
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
    return((3*n0[0]*n0[1]*n0[2] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 9*bs));
}

/*
    Full Multigrid solver.  See Numerical Recipes (second edition) for more
    information
*/
void fmg3_noa(int n0[], float *b0, int rtype, double param0[], int c, int nit,
          float *u0, float *scratch)
{
    int i, j, ng, bs;
     int n[32][3], m[32];
    float *bo[32], *b[32], *u[32], *res, *rbuf;
    double param[32][6];
    void (*relax)(), (*LtLf)();
    double (*sumsq)();

    if (rtype == 0)
    {
        relax = relax_le_noa;
        LtLf  = LtLf_le;
        sumsq = sumsq_le_noa;
    }
    else
        return;

    bo[0]   = b0;
    b[0]    = b0;
    u[0]    = u0;
    n[0][0] = n0[0];
    n[0][1] = n0[1];
    n[0][2] = n0[2];
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
        m[j]    = n[j][0]*n[j][1]*n[j][2];
        ng ++;
        bs += m[j];
        if ((n[j][0]<2) && (n[j][1]<2) && (n[j][2]<2))
            break;
    }

    res    = scratch;
    rbuf   = scratch + 3*m[0];
    bo[1]  = scratch + 3*m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1];
    b[1]   = scratch + 3*m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 3*bs;
    u[1]   = scratch + 3*m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 6*bs;

    for(j=2; j<ng; j++)
    {
        bo[j] = bo[j-1]+3*m[j-1];
        b[j]  =  b[j-1]+3*m[j-1];
        u[j]  =  u[j-1]+3*m[j-1];
    }

    for(j=1; j<ng; j++)
    {
        restrict(3,n[j-1],bo[j-1],n[j],bo[j],rbuf);

        param[j][0] = param0[0]*(double)n[j][0]/n0[0];
        param[j][1] = param0[1]*(double)n[j][1]/n0[1];
        param[j][2] = param0[2]*(double)n[j][2]/n0[2];
        param[j][3] = param[0][3];
        param[j][4] = param[0][4];
        param[j][5] = param[0][5];
    }

    if (u[0][0]==0)
    {
     /* solve33_noa(bo[ng-1], param0[5], u[ng-1]); */
        relax(n[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]);

        for(j=ng-2; j>=0; j--)
        {
            int jc;
            prolong(3,n[j+1],u[j+1],n[j],u[j],rbuf);
            if(j>0) copy(3*m[j],bo[j],b[j]);
            for(jc=0; jc<c; jc++)
            {
                int jj;
                for(jj=j; jj<ng-1; jj++)
                {
                    relax(n[jj], b[jj], param[jj], nit, u[jj]);
                    LtLf(n[jj], u[jj], param[jj], res);
                    for(i=0; i<3*m[jj]; i++)
                        res[i] = b[jj][i] - res[i];
                    restrict(3,n[jj],res,n[jj+1],b[jj+1],rbuf);
                    zeros(3*m[jj+1],u[jj+1]);
                }
             /* solve33_noa(b[ng-1], param0[5], u[ng-1]); */
                relax(n[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]);
                for(jj=ng-2; jj>=j; jj--)
                {
                    prolong(3,n[jj+1],u[jj+1],n[jj],res,rbuf);
                    addto(3*m[jj], u[jj], res);
                    relax(n[jj], b[jj], param[jj], nit, u[jj]);
                }
            }
        }
    }
    else
    {
        int jc;
        for(j=1; j<ng; j++)
            restrict(3,n[j-1],u[j-1],n[j],u[j],rbuf);

        for(jc=0; jc<c; jc++)
        {
            int jj;
            for(jj=0; jj<ng-1; jj++)
            {
                relax(n[jj], b[jj], param[jj], nit, u[jj]);
                LtLf(n[jj], u[jj], param[jj], res);
                for(i=0; i<3*m[jj]; i++)
                    res[i] = b[jj][i] - res[i];
                restrict(3,n[jj],res,n[jj+1],b[jj+1],rbuf);
                zeros(3*m[jj+1],u[jj+1]);
            }
         /* solve33_noa(b[ng-1], param0[5], u[ng-1]); */
            relax(n[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]);
            for(jj=ng-2; jj>=0; jj--)
            {
                prolong(3,n[jj+1],u[jj+1],n[jj],res,rbuf);
                addto(3*m[jj], u[jj], res);
                relax(n[jj], b[jj], param[jj], nit, u[jj]);
            }
        }
    }
}

