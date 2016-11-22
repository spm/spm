/* $Id: shoot_diffeo3d.c 6945 2016-11-22 10:29:18Z john $ */
/* (c) John Ashburner (2011) */

#include <mex.h>
#include <math.h>
#include <stdio.h>
#include "shoot_optim3d.h"
#include "shoot_expm3.h"
#include "shoot_boundary.h"

extern double   log(double x);
extern double   exp(double x);
#define LOG(x) (((x)>0) ? log(x+0.001): -6.9078)

/*
 * Lie Bracket
 * C = [A,B]
 */
void bracket(mwSize dm[], float *A, float *B, float *C)
{
    float *Ax, *Ay, *Az;
    float *Bx, *By, *Bz;
    float *Cx, *Cy, *Cz;
    mwSize i, j, k, mm = dm[0]*dm[1]*dm[2];

    Ax = A;
    Ay = A + mm;
    Az = A + mm*2;

    Bx = B;
    By = B + mm;
    Bz = B + mm*2;

    Cx = C;
    Cy = C + mm;
    Cz = C + mm*2;

    for(k=0; k<dm[2]; k++)
    {
        for(j=0; j<dm[1]; j++)
        {
            mwSize o1, oi1, opj1, omj1, opk1, omk1;
            o1   = dm[0]*(j+dm[1]*k);
            oi1  = dm[0]*(j+dm[1]*k);
            opj1 = dm[0]*(bound(j+1,dm[1])+dm[1]*k);
            omj1 = dm[0]*(bound(j-1,dm[1])+dm[1]*k);
            opk1 = dm[0]*(j+dm[1]*bound(k+1,dm[2]));
            omk1 = dm[0]*(j+dm[1]*bound(k-1,dm[2]));

            for(i=0; i<dm[0]; i++)
            {
                mwSize o, opi, omi, opj, omj, opk, omk;
                double j00, j01, j02,  j10, j11, j12,  j20, j21, j22;
                double tx, ty, tz,  cx1, cy1, cz1,  cx2, cy2, cz2;

                o   = i+o1;
                opi = bound(i+1,dm[0])+oi1;
                omi = bound(i-1,dm[0])+oi1;
                opj = i+opj1;
                omj = i+omj1;
                opk = i+opk1;
                omk = i+omk1;

                tx = Ax[o];
                ty = Ay[o];
                tz = Az[o];

                j00 = (Bx[opi]-Bx[omi])/2.0;
                j01 = (By[opi]-By[omi])/2.0;
                j02 = (Bz[opi]-Bz[omi])/2.0;

                j10 = (Bx[opj]-Bx[omj])/2.0;
                j11 = (By[opj]-By[omj])/2.0;
                j12 = (Bz[opj]-Bz[omj])/2.0;

                j20 = (Bx[opk]-Bx[omk])/2.0;
                j21 = (By[opk]-By[omk])/2.0;
                j22 = (Bz[opk]-Bz[omk])/2.0;

                cx1 = tx*j00+ty*j10+tz*j20;
                cy1 = tx*j01+ty*j11+tz*j21;
                cz1 = tx*j02+ty*j12+tz*j22;

                tx = Bx[o];
                ty = By[o];
                tz = Bz[o];

                j00 = (Ax[opi]-Ax[omi])/2.0;
                j01 = (Ay[opi]-Ay[omi])/2.0;
                j02 = (Az[opi]-Az[omi])/2.0;

                j10 = (Ax[opj]-Ax[omj])/2.0;
                j11 = (Ay[opj]-Ay[omj])/2.0;
                j12 = (Az[opj]-Az[omj])/2.0;

                j20 = (Ax[opk]-Ax[omk])/2.0;
                j21 = (Ay[opk]-Ay[omk])/2.0;
                j22 = (Az[opk]-Az[omk])/2.0;

                cx2 = tx*j00+ty*j10+tz*j20;
                cy2 = tx*j01+ty*j11+tz*j21;
                cz2 = tx*j02+ty*j12+tz*j22;

                Cx[o] = cx2-cx1;
                Cy[o] = cy2-cy1;
                Cz[o] = cz2-cz1;
            }
        }
    }
}

/*
 * Composition operations, possibly along with Jacobian matrices
 */
static void composition_stuff(mwSize dm[], mwSize mm,
                              float *B, float *JB, float *A, float *JA,
                              float *C, float *JC, int flag)
{
    float *Ax, *Ay, *Az, *JA00, *JA01, *JA02,  *JA10, *JA11, *JA12,  *JA20, *JA21, *JA22;
    float *Bx, *By, *Bz, *JB00, *JB01, *JB02,  *JB10, *JB11, *JB12,  *JB20, *JB21, *JB22;
    float *Cx, *Cy, *Cz, *JC00, *JC01, *JC02,  *JC10, *JC11, *JC12,  *JC20, *JC21, *JC22;
    mwSize i, mmb = dm[0]*dm[1]*dm[2];

    /* Does not yet work properly if dimensions of A and B are not identical.
       Still need to figure out why not. */

    Ax   =  A;
    Ay   =  A+mm;
    Az   =  A+mm*2;
    Bx   =  B;
    By   =  B+mmb;
    Bz   =  B+mmb*2;
    Cx   =  C;
    Cy   =  C+mm;
    Cz   =  C+mm*2;

    if (JC && !flag)
    {
        JA00 = JA+mm*0; JA01 = JA+mm*1; JA02 = JA+mm*2;
        JA10 = JA+mm*3; JA11 = JA+mm*4; JA12 = JA+mm*5;
        JA20 = JA+mm*6; JA21 = JA+mm*7; JA22 = JA+mm*8;

        JB00 = JB+mmb*0; JB01 = JB+mmb*1; JB02 = JB+mmb*2;
        JB10 = JB+mmb*3; JB11 = JB+mmb*4; JB12 = JB+mmb*5;
        JB20 = JB+mmb*6; JB21 = JB+mmb*7; JB22 = JB+mmb*8;

        JC00 = JC+mm*0; JC01 = JC+mm*1; JC02 = JC+mm*2;
        JC10 = JC+mm*3; JC11 = JC+mm*4; JC12 = JC+mm*5;
        JC20 = JC+mm*6; JC21 = JC+mm*7; JC22 = JC+mm*8;
    }

    for(i=0; i<mm; i++)
    {
        double x, y, z;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
        mwSize o000, o100, o010, o110, o001, o101, o011, o111;
        mwSize tmpz, tmpy, n;

        x    = Ax[i]-1.0;
        y    = Ay[i]-1.0;
        z    = Az[i]-1.0;
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
        o000  = ix +tmpy;
        o100  = ix1+tmpy;
        tmpy  = dm[0]*(iy1 + tmpz);
        o010  = ix +tmpy;
        o110  = ix1+tmpy;
        tmpz  = dm[1]*iz1;
        tmpy  = dm[0]*(iy + tmpz);
        o001  = ix +tmpy;
        o101  = ix1+tmpy;
        tmpy  = dm[0]*(iy1 + tmpz);
        o011  = ix +tmpy;
        o111  = ix1+tmpy;

        k000  = Bx[o000]-1.0;
        k100  = Bx[o100]-1.0;
        k010  = Bx[o010]-1.0;
        k110  = Bx[o110]-1.0;
        k001  = Bx[o001]-1.0;
        k101  = Bx[o101]-1.0;
        k011  = Bx[o011]-1.0;
        k111  = Bx[o111]-1.0;

        n     = dm[0];
        k100 -= floor((k100-k000)/n+0.5)*n;
        k010 -= floor((k010-k000)/n+0.5)*n;
        k110 -= floor((k110-k000)/n+0.5)*n;
        k001 -= floor((k001-k000)/n+0.5)*n;
        k101 -= floor((k101-k000)/n+0.5)*n;
        k011 -= floor((k011-k000)/n+0.5)*n;
        k111 -= floor((k111-k000)/n+0.5)*n;
        Cx[i] = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
              + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 + 1.0;

        k000  = By[o000]-1.0;
        k100  = By[o100]-1.0;
        k010  = By[o010]-1.0;
        k110  = By[o110]-1.0;
        k001  = By[o001]-1.0;
        k101  = By[o101]-1.0;
        k011  = By[o011]-1.0;
        k111  = By[o111]-1.0;

        n     = dm[1];
        k100 -= floor((k100-k000)/n+0.5)*n;
        k010 -= floor((k010-k000)/n+0.5)*n;
        k110 -= floor((k110-k000)/n+0.5)*n;
        k001 -= floor((k001-k000)/n+0.5)*n;
        k101 -= floor((k101-k000)/n+0.5)*n;
        k011 -= floor((k011-k000)/n+0.5)*n;
        k111 -= floor((k111-k000)/n+0.5)*n;
        Cy[i] = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
              + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 + 1.0;

        k000  = Bz[o000]-1.0;
        k100  = Bz[o100]-1.0;
        k010  = Bz[o010]-1.0;
        k110  = Bz[o110]-1.0;
        k001  = Bz[o001]-1.0;
        k101  = Bz[o101]-1.0;
        k011  = Bz[o011]-1.0;
        k111  = Bz[o111]-1.0;

        n     = dm[2];
        k100 -= floor((k100-k000)/n+0.5)*n;
        k010 -= floor((k010-k000)/n+0.5)*n;
        k110 -= floor((k110-k000)/n+0.5)*n;
        k001 -= floor((k001-k000)/n+0.5)*n;
        k101 -= floor((k101-k000)/n+0.5)*n;
        k011 -= floor((k011-k000)/n+0.5)*n;
        k111 -= floor((k111-k000)/n+0.5)*n;
        Cz[i] = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
              + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 + 1.0;

        if (JC)
        {
            if (flag==0)
            {
                float *ptr;
                double ja0, ja1, ja2;
                double jb[3][3];

                ptr      = JB00;
                jb[0][0] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;
                ptr      = JB10;
                jb[1][0] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;
                ptr      = JB20;
                jb[2][0] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;

                ptr      = JB01;
                jb[0][1] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;
                ptr      = JB11;
                jb[1][1] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;
                ptr      = JB21;
                jb[2][1] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;

                ptr      = JB02;
                jb[0][2] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;

                ptr      = JB12;
                jb[1][2] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;

                ptr      = JB22;
                jb[2][2] = ((ptr[o000]*dx2 + ptr[o100]*dx1)*dy2 + (ptr[o010]*dx2 + ptr[o110]*dx1)*dy1)*dz2
                         + ((ptr[o001]*dx2 + ptr[o101]*dx1)*dy2 + (ptr[o011]*dx2 + ptr[o111]*dx1)*dy1)*dz1;

                ja0     = JA00[i];
                ja1     = JA01[i];
                ja2     = JA02[i];
                JC00[i] = jb[0][0]*ja0 + jb[1][0]*ja1 + jb[2][0]*ja2;
                JC01[i] = jb[0][1]*ja0 + jb[1][1]*ja1 + jb[2][1]*ja2;
                JC02[i] = jb[0][2]*ja0 + jb[1][2]*ja1 + jb[2][2]*ja2;

                ja0     = JA10[i];
                ja1     = JA11[i];
                ja2     = JA12[i];
                JC10[i] = jb[0][0]*ja0 + jb[1][0]*ja1 + jb[2][0]*ja2;
                JC11[i] = jb[0][1]*ja0 + jb[1][1]*ja1 + jb[2][1]*ja2;
                JC12[i] = jb[0][2]*ja0 + jb[1][2]*ja1 + jb[2][2]*ja2;

                ja0     = JA20[i];
                ja1     = JA21[i];
                ja2     = JA22[i];
                JC20[i] = jb[0][0]*ja0 + jb[1][0]*ja1 + jb[2][0]*ja2;
                JC21[i] = jb[0][1]*ja0 + jb[1][1]*ja1 + jb[2][1]*ja2;
                JC22[i] = jb[0][2]*ja0 + jb[1][2]*ja1 + jb[2][2]*ja2;
            }
            else
            {
                double jb;
                jb    = ((JB[o000]*dx2 + JB[o100]*dx1)*dy2 + (JB[o010]*dx2 + JB[o110]*dx1)*dy1)*dz2
                      + ((JB[o001]*dx2 + JB[o101]*dx1)*dy2 + (JB[o011]*dx2 + JB[o111]*dx1)*dy1)*dz1;
                JC[i] = jb * JA[i];
            }
        }
    }
}

/*
 * Composition operation
 * C(Id) = B(A(Id))
 */
void composition(mwSize dm[], mwSize mm, float *B, float *A, float *C)
{
    composition_stuff(dm, mm, B, 0, A, 0, C, 0,0);
}

/*
 * Composition operation, along with Jacobian matrices
 * C(Id)  =  B(A(Id))
 * JC(Id) = JB(A(Id))*JA(Id) ?
 */
void composition_jacobian(mwSize dm[], mwSize mm, float *B, float *JB, float *A, float *JA, float *C, float *JC)
{
    composition_stuff(dm, mm, B, JB, A, JA, C, JC, 0);
}

/*
 * Composition operation, along with Jacobian determinants
 * C(Id)  =  B(A(Id))
 * JC(Id) = JB(A(Id))*JA(Id)
 */
void composition_jacdet(mwSize dm[], mwSize mm, float *B, float *JB, float *A, float *JA, float *C, float *JC)
{
    composition_stuff(dm, mm, B, JB, A, JA, C, JC, 1);
}

static void def2jac_wrap(mwSize dm[], float *Y, float *J, mwSignedIndex s, int code)
{
    mwSignedIndex i, j, k, k0, k2, mm;

    if (s>=0 && s<dm[2])
    {
        k0 = s;
        k2 = s+1;
    }
    else
    {
        k0 = 0;
        k2 = dm[2];
    }

    mm = dm[0]*dm[1]*(k2-k0);

    for(k=k0; k<k2; k++)
    {
        mwSignedIndex km, kp;
        km = (bound(k-1,dm[2])-k)*dm[1]*dm[0];
        kp = (bound(k+1,dm[2])-k)*dm[1]*dm[0];

        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex jm, jp;
            float *y1, *y2, *y3, *dp;
            jm = (bound(j-1,dm[1])-j)*dm[0];
            jp = (bound(j+1,dm[1])-j)*dm[0];
            y1 = Y + dm[0]*(j+dm[1]*k);
            y2 = Y + dm[0]*(j+dm[1]*(k+dm[2]));
            y3 = Y + dm[0]*(j+dm[1]*(k+dm[2]*2));
            dp = J + dm[0]*(j+dm[1]*k);

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im, ip;
                double j11, j12, j13, j21, j22, j23, j31, j32, j33;
                im = bound(i-1,dm[0])-i;
                ip = bound(i+1,dm[0])-i;

                j11 = (y1[i+ip]-y1[i+im])/dm[0]; j11 = 0.5*(j11 - floor(j11+0.5))*dm[0];
                j21 = (y2[i+ip]-y2[i+im])/dm[0]; j21 = 0.5*(j21 - floor(j21+0.5))*dm[0];
                j31 = (y3[i+ip]-y3[i+im])/dm[0]; j31 = 0.5*(j31 - floor(j31+0.5))*dm[0];

                j12 = (y1[i+jp]-y1[i+jm])/dm[1]; j12 = 0.5*(j12 - floor(j12+0.5))*dm[1];
                j22 = (y2[i+jp]-y2[i+jm])/dm[1]; j22 = 0.5*(j22 - floor(j22+0.5))*dm[1];
                j32 = (y3[i+jp]-y3[i+jm])/dm[1]; j32 = 0.5*(j32 - floor(j32+0.5))*dm[1];

                if (dm[2]>1)
                {
                    j13 = (y1[i+kp]-y1[i+km])/dm[2]; j13 = 0.5*(j13 - floor(j13+0.5))*dm[2];
                    j23 = (y2[i+kp]-y2[i+km])/dm[2]; j23 = 0.5*(j23 - floor(j23+0.5))*dm[2];
                    j33 = (y3[i+kp]-y3[i+km])/dm[2]; j33 = 0.5*(j33 - floor(j33+0.5))*dm[2];
                }
                else
                {
                    j13 = 0.0;
                    j23 = 0.0;
                    j33 = 1.0;
                }

                if (!code)
                {
                    dp[i] = j11*(j22*j33-j23*j32)+j21*(j13*j32-j12*j33)+j31*(j12*j23-j13*j22);
                }
                else
                {
                    dp[i     ] = j11; dp[i+mm  ] = j21; dp[i+mm*2] = j31;
                    dp[i+mm*3] = j12; dp[i+mm*4] = j22; dp[i+mm*5] = j32;
                    dp[i+mm*6] = j13; dp[i+mm*7] = j23; dp[i+mm*8] = j33;
                }
            }
        }
    }
}

static void def2jac_neuman(mwSize dm[], float *Y, float *J, mwSignedIndex s, int code)
{
    mwSignedIndex k, k0, k2, mm;

    if (s>=0 && s<dm[2])
    {
        k0 = s;
        k2 = s;
    }
    else
    {
        k0 = 0;
        k2 = dm[2];
    }
    mm = dm[0]*dm[1]*(k2-k0);

    for(k=k0; k<k2; k++)
    {
        mwSignedIndex j;
        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex i;
            float *y1, *y2, *y3, *dp;
            y1 = Y + dm[0]*(j+dm[1]*k);
            y2 = Y + dm[0]*(j+dm[1]*(k+dm[2]));
            y3 = Y + dm[0]*(j+dm[1]*(k+dm[2]*2));
            dp = J + dm[0]*(j+dm[1]*k);

            for(i=0; i<dm[0]; i++)
            {
                double j11, j12, j13, j21, j22, j23, j31, j32, j33;

                if ((i==0) || (i==dm[0]-1))
                {
                    j11 = 1.0;
                    j21 = 0.0;
                    j31 = 0.0;
                }
                else
                {
                    j11 = (y1[i+1]-y1[i-1])*0.5;
                    j21 = (y2[i+1]-y2[i-1])*0.5;
                    j31 = (y3[i+1]-y3[i-1])*0.5;
                }

                if ((j==0) || (j==dm[1]-1))
                {
                    j12 = 0.0;
                    j22 = 1.0;
                    j32 = 0.0;
                }
                else
                {
                    mwSignedIndex op = i+dm[0], om = i-dm[0];
                    j12 = (y1[op]-y1[om])*0.5;
                    j22 = (y2[op]-y2[om])*0.5;
                    j32 = (y3[op]-y3[om])*0.5;
                }

                if ((k==0) || (k==dm[2]-1))
                {
                    j13 = 0.0;
                    j23 = 0.0;
                    j33 = 1.0;
                }
                else
                {
                    mwSignedIndex op = i+dm[0]*dm[1], om = i-dm[0]*dm[1];
                    j13 = (y1[op]-y1[om])*0.5; 
                    j23 = (y2[op]-y2[om])*0.5;
                    j33 = (y3[op]-y3[om])*0.5;
                }

                if (!code)
                {
                    dp[i] = j11*(j22*j33-j23*j32)+j21*(j13*j32-j12*j33)+j31*(j12*j23-j13*j22);
                }
                else
                {
                    dp[i     ] = j11; dp[i+mm  ] = j21; dp[i+mm*2] = j31;
                    dp[i+mm*3] = j12; dp[i+mm*4] = j22; dp[i+mm*5] = j32;
                    dp[i+mm*6] = j13; dp[i+mm*7] = j23; dp[i+mm*8] = j33;
                }
            }
        }
    }
}

void def2det(mwSize dm[], float *Y, float *J, mwSignedIndex s)
{
    if (get_bound())
        def2jac_neuman(dm, Y, J, s, 0);
    else
        def2jac_wrap(dm, Y, J, s, 0);
}

void def2jac(mwSize dm[], float *Y, float *J, mwSignedIndex s)
{
    if (get_bound())
        def2jac_neuman(dm, Y, J, s, 1);
    else
        def2jac_wrap(dm, Y, J, s, 1);
}

/*
 * Sample a point
 * s = f(x,y,z)
 */
double samp(mwSize dm[], float f[], double x, double y, double z)
{
    mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
    mwSize o000, o100, o010, o110, o001, o101, o011, o111;
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
    o000  = ix +tmpy;
    o100  = ix1+tmpy;
    tmpy  = dm[0]*(iy1 + tmpz);
    o010  = ix +tmpy;
    o110  = ix1+tmpy;
    tmpz  = dm[1]*iz1;
    tmpy  = dm[0]*(iy + tmpz);
    o001  = ix +tmpy;
    o101  = ix1+tmpy;
    tmpy  = dm[0]*(iy1 + tmpz);
    o011  = ix +tmpy;
    o111  = ix1+tmpy;

    return( ((f[o000]*dx2 + f[o100]*dx1)*dy2 + (f[o010]*dx2 + f[o110]*dx1)*dy1)*dz2
          + ((f[o001]*dx2 + f[o101]*dx1)*dy2 + (f[o011]*dx2 + f[o111]*dx1)*dy1)*dz1 );
}

/* Sample n points
 * s1 = f1(x,y,z)
 * s2 = f2(x,y,z)
 */
void sampn(mwSize dm[], float f[], mwSize n, mwSize mm, double x, double y, double z, double v[])
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
    o000  = ix +tmpy;
    o100  = ix1+tmpy;
    tmpy  = dm[0]*(iy1 + tmpz);
    o010  = ix +tmpy;
    o110  = ix1+tmpy;
    tmpz  = dm[1]*iz1;
    tmpy  = dm[0]*(iy + tmpz);
    o001  = ix +tmpy;
    o101  = ix1+tmpy;
    tmpy  = dm[0]*(iy1 + tmpz);
    o011  = ix +tmpy;
    o111  = ix1+tmpy;

    for(j=0; j<n; j++, f += mm)
    {
        v[j] = ((f[o000]*dx2 + f[o100]*dx1)*dy2 + (f[o010]*dx2 + f[o110]*dx1)*dy1)*dz2
             + ((f[o001]*dx2 + f[o101]*dx1)*dy2 + (f[o011]*dx2 + f[o111]*dx1)*dy1)*dz1;
    }
}

/* Rather than sample from an image according to a deformation,
 * it is also possible to push voxels from one image into
 * another according to the inverse of the deformation.
 * Note that the result is a noisy version of a Jacobian "modulated"
 * image.
 */
void push(mwSize dm[], mwSize m, mwSize n, float def[], float pf[], float po[], float so[])
{
    mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
    mwSize   i, j, mm, tmpz, tmpy;
    float *px, *py, *pz;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    px = def;
    py = def+m;
    pz = def+m*2;
    mm = dm[0]*dm[1]*dm[2];

    for(i=0; i<m; i++)
    {
        double x, y, z;

        if (mxIsFinite(pf[i]) && mxIsFinite(px[i]) && mxIsFinite(py[i]) && mxIsFinite(pz[i]))
        {
            x    = px[i]-1.0; /* Subtract 1 because of MATLAB indexing */
            y    = py[i]-1.0;
            z    = pz[i]-1.0;

            /* Check range and avoid inserting values outside the FOV. */
            if (x>=0 && x<dm[0]-1 && y>=0 && y<dm[1]-1 && z>=0 && z<dm[2]-1)
            {
                /* A faster function fo voxels that are safely inside the FOV */
                mwSize o000, o100, o010, o110, o001, o101, o011, o111;
                float w000, w100, w010, w110, w001, w101, w011, w111;
                ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
                iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
                iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;

                /* Weights for trilinear interpolation */
                w000 = dx2*dy2*dz2;
                w100 = dx1*dy2*dz2;
                w010 = dx2*dy1*dz2;
                w110 = dx1*dy1*dz2;
                w001 = dx2*dy2*dz1;
                w101 = dx1*dy2*dz1;
                w011 = dx2*dy1*dz1;
                w111 = dx1*dy1*dz1;

                ix1  = ix+1;
                iy1  = iy+1;
                iz1  = iz+1;

                /* Neighbouring voxels used for trilinear interpolation */
                tmpz  = dm[1]*iz;
                tmpy  = dm[0]*(iy + tmpz);
                o000  = ix +tmpy;
                o100  = ix1+tmpy;
                tmpy  = dm[0]*(iy1 + tmpz);
                o010  = ix +tmpy;
                o110  = ix1+tmpy;
                tmpz  = dm[1]*iz1;
                tmpy  = dm[0]*(iy + tmpz);
                o001  = ix +tmpy;
                o101  = ix1+tmpy;
                tmpy  = dm[0]*(iy1 + tmpz);
                o011  = ix +tmpy;
                o111  = ix1+tmpy;

                for (j=0; j<n; j++)
                {
                    /* Increment the images themselves */
                    float *pj = po+mm*j;
                    float  f  = pf[i+j*m];
                    pj[o000] += f*w000;
                    pj[o100] += f*w100;
                    pj[o010] += f*w010;
                    pj[o110] += f*w110;
                    pj[o001] += f*w001;
                    pj[o101] += f*w101;
                    pj[o011] += f*w011;
                    pj[o111] += f*w111;
                }

                if (so!=(float *)0)
                {
                    /* Increment an image containing the number of voxels added */
                    so[o000] += w000;
                    so[o100] += w100;
                    so[o010] += w010;
                    so[o110] += w110;
                    so[o001] += w001;
                    so[o101] += w101;
                    so[o011] += w011;
                    so[o111] += w111;
                }
            }
            else if ((x>=-1) && (x<dm[0]) && (y>=-1) && (y<dm[1]) && (z>=-1) && (z<dm[2]))
            {
                /* A slower function for voxels at the edge of the field of view */
                mwSize o[8], nn=0, k;
                float w[8];

                ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
                iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
                iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;
                ix1  = ix+1;
                iy1  = iy+1;
                iz1  = iz+1;
                if (iz>=0)
                {
                    tmpz  = dm[1]*iz;
                    if (iy>=0)
                    {
                        tmpy  = dm[0]*(iy + tmpz);
                        if (ix>=0)
                        {
                            o[nn] = ix+tmpy;
                            w[nn] = dx2*dy2*dz2;
                            nn++;
                        }
                        if (ix1<dm[0])
                        {
                            o[nn] = ix1+tmpy;
                            w[nn] = dx1*dy2*dz2;
                            nn++;
                        }
                    }
                    if (iy1<dm[1])
                    {
                        tmpy  = dm[0]*(iy1 + tmpz);
                        if (ix>=0)
                        {
                            o[nn] = ix+tmpy;
                            w[nn] = dx2*dy1*dz2;
                            nn++;
                        }
                        if (ix1<dm[0])
                        {
                            o[nn] = ix1+tmpy;
                            w[nn] = dx1*dy1*dz2;
                            nn++;
                        }
                    }
                }
                if (iz1<dm[2])
                {
                    tmpz  = dm[1]*iz1;
                    if (iy>=0)
                    {
                        tmpy  = dm[0]*(iy + tmpz);
                        if (ix>=0)
                        {
                            o[nn] = ix +tmpy;
                            w[nn] = dx2*dy2*dz1;
                            nn++;
                        }
                        if (ix1<dm[0])
                        {
                            o[nn] = ix1+tmpy;
                            w[nn] = dx1*dy2*dz1;
                            nn++;
                        }
                    }
                    if (iy1<dm[1])
                    {
                        tmpy  = dm[0]*(iy1 + tmpz);
                        if (ix>=0)
                        {
                            o[nn] = ix +tmpy;
                            w[nn] = dx2*dy1*dz1;
                            nn++;
                        }
                        if (ix1<dm[0])
                        {
                            o[nn] = ix1+tmpy;
                            w[nn] = dx1*dy1*dz1;
                            nn++;
                        }
                    }
                }
                if (so!=(float *)0)
                {
                    for(k=0; k<nn; k++)
                        so[o[k]] += w[k];
                }

                for (j=0; j<n; j++)
                {
                    float *pj = po+mm*j;
                    float  f  = pf[i+j*m];
                    for(k=0; k<nn; k++)
                        pj[o[k]] += f*w[k];
                }
            }
        }
    }
}

/* Same as above, but with circulant boundary conditions */
void pushc(mwSize dm[], mwSize m, mwSize n, float def[], float pf[], float po[], float so[])
{
    mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
    mwSize i, j, mm, tmpz, tmpy;
    float *px, *py, *pz;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    px = def;
    py = def+m;
    pz = def+m*2;
    mm = dm[0]*dm[1]*dm[2];

    for(i=0; i<m; i++)
    {
        double x, y, z;

        if (mxIsFinite(pf[i]) && mxIsFinite(px[i]) && mxIsFinite(py[i]) && mxIsFinite(pz[i]))
        {
            mwSize o000, o100, o010, o110, o001, o101, o011, o111;
            float w000, w100, w010, w110, w001, w101, w011, w111;

            x    = px[i]-1.0; /* Subtract 1 because of MATLAB indexing */
            y    = py[i]-1.0;
            z    = pz[i]-1.0;

            ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
            iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
            iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;

            /* Weights for trilinear interpolation */
            w000 = dx2*dy2*dz2;
            w100 = dx1*dy2*dz2;
            w010 = dx2*dy1*dz2;
            w110 = dx1*dy1*dz2;
            w001 = dx2*dy2*dz1;
            w101 = dx1*dy2*dz1;
            w011 = dx2*dy1*dz1;
            w111 = dx1*dy1*dz1;

            ix   = bound(ix, dm[0]);
            iy   = bound(iy, dm[1]);
            iz   = bound(iz, dm[2]);
            ix1  = bound(ix+1, dm[0]);
            iy1  = bound(iy+1, dm[1]);
            iz1  = bound(iz+1, dm[2]);

            /* Neighbouring voxels used for trilinear interpolation */
            tmpz  = dm[1]*iz;
            tmpy  = dm[0]*(iy + tmpz);
            o000  = ix +tmpy;
            o100  = ix1+tmpy;
            tmpy  = dm[0]*(iy1 + tmpz);
            o010  = ix +tmpy;
            o110  = ix1+tmpy;
            tmpz  = dm[1]*iz1;
            tmpy  = dm[0]*(iy + tmpz);
            o001  = ix +tmpy;
            o101  = ix1+tmpy;
            tmpy  = dm[0]*(iy1 + tmpz);
            o011  = ix +tmpy;
            o111  = ix1+tmpy;

            for (j=0; j<n; j++)
            {
                /* Increment the images themselves */
                float *pj = po+mm*j;
                float  f  = pf[i+j*m];
                pj[o000] += f*w000;
                pj[o100] += f*w100;
                pj[o010] += f*w010;
                pj[o110] += f*w110;
                pj[o001] += f*w001;
                pj[o101] += f*w101;
                pj[o011] += f*w011;
                pj[o111] += f*w111;
            }

            if (so!=(float *)0)
            {
                /* Increment an image containing the number of voxels added */
                so[o000] += w000;
                so[o100] += w100;
                so[o010] += w010;
                so[o110] += w110;
                so[o001] += w001;
                so[o101] += w101;
                so[o011] += w011;
                so[o111] += w111;
            }
        }
    }
}

/* Similar to above, except with a multiplication by the inverse of the Jacobians.
   This is used for geodesic shooting */
void pushc_grads(mwSize dmo[], mwSize dm[], float def[], float J[], float pf[], float po[])
{
    mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
    mwSize i2, i, mo, my;
    float *px, *py, *pz;
    float *pj11, *pj12, *pj13, *pj21, *pj22, *pj23, *pj31, *pj32, *pj33;

    my = dm[0]*dm[1]*dm[2];
    px = def;
    py = def+my;
    pz = def+my*2;
    mo = dmo[0]*dmo[1]*dmo[2];

    if (J!=(float *)0)
    {
        pj11 = J;
        pj21 = J+my;
        pj31 = J+my*2;
        pj12 = J+my*3;
        pj22 = J+my*4;
        pj32 = J+my*5;
        pj13 = J+my*6;
        pj23 = J+my*7;
        pj33 = J+my*8;
    }

    for(i2=0, i=0; i2<dm[2]; i2++)
    {
        mwSize i1;
        mwSignedIndex i2m, i2p;

        i2m = (bound(i2-1,dm[2])-i2)*dm[1]*dm[0];
        i2p = (bound(i2+1,dm[2])-i2)*dm[1]*dm[0];

        for(i1=0; i1<dm[1]; i1++)
        {
            mwSize i0;
            mwSignedIndex i1m, i1p;

            i1m = (bound(i1-1,dm[1])-i1)*dm[0];
            i1p = (bound(i1+1,dm[1])-i1)*dm[0];

            for(i0=0; i0<dm[0]; i0++, i++, px++, py++, pz++)
            {
                if (mxIsFinite(pf[i]) && mxIsFinite(*px) && mxIsFinite(*py) && mxIsFinite(*pz))
                {
                    mwSize j, tmpz, tmpy;
                    mwSize o000, o100, o010, o110, o001, o101, o011, o111;
                    float  w000, w100, w010, w110, w001, w101, w011, w111;
                    float j11, j12, j13, j21, j22, j23, j31, j32, j33;
                    float ij11, ij12, ij13, ij21, ij22, ij23, ij31, ij32, ij33, dj;
                    float f1, f2, f3;
                    double x, y, z;
                    double rf[3];
                    double dx1, dx2, dy1, dy2, dz1, dz2;

                    x    = *px-1.0; /* Subtract 1 because of MATLAB indexing */
                    y    = *py-1.0;
                    z    = *pz-1.0;

                    if (J!=(float *)0)
                    {
                        /* If Jacobians are passed, use them */
                        j11  = pj11[i]; j12  = pj12[i]; j13  = pj13[i];
                        j21  = pj21[i]; j22  = pj22[i]; j23  = pj23[i];
                        j31  = pj31[i]; j32  = pj32[i]; j33  = pj33[i];
                    }
                    else
                    {
                        /* Otherwise compute Jacobians from deformation */
                        if (get_bound()==0)
                        {
                            /* Circulant boundary */
                            mwSignedIndex i0m, i0p;
                            i0m = bound(i0-1,dm[0])-i0;
                            i0p = bound(i0+1,dm[0])-i0;

                            j11 = (px[i0p]-px[i0m])/dm[0]; j11 = 0.5*(j11 - floor(j11+0.5))*dm[0];
                            j21 = (py[i0p]-py[i0m])/dm[0]; j21 = 0.5*(j21 - floor(j21+0.5))*dm[0];
                            j31 = (pz[i0p]-pz[i0m])/dm[0]; j31 = 0.5*(j31 - floor(j31+0.5))*dm[0];

                            j12 = (px[i1p]-px[i1m])/dm[1]; j12 = 0.5*(j12 - floor(j12+0.5))*dm[1];
                            j22 = (py[i1p]-py[i1m])/dm[1]; j22 = 0.5*(j22 - floor(j22+0.5))*dm[1];
                            j32 = (pz[i1p]-pz[i1m])/dm[1]; j32 = 0.5*(j32 - floor(j32+0.5))*dm[1];

                            if (dm[2]>1)
                            {
                                j13 = (px[i2p]-px[i2m])/dm[2]; j13 = 0.5*(j13 - floor(j13+0.5))*dm[2];
                                j23 = (py[i2p]-py[i2m])/dm[2]; j23 = 0.5*(j23 - floor(j23+0.5))*dm[2];
                                j33 = (pz[i2p]-pz[i2m])/dm[2]; j33 = 0.5*(j33 - floor(j33+0.5))*dm[2];
                            }
                            else
                            {
                                j13 = 0.0;
                                j23 = 0.0;
                                j33 = 1.0;
                            }
                        }
                        else
                        {
                            /* Neumann boundary - only discontinuities at edges */
                            if ((i0==0) || (i0==dm[0]-1))
                            {
                                j11 = 1.0;
                                j21 = 0.0;
                                j31 = 0.0;
                            }
                            else
                            {
                                j11 = (px[1]-px[-1])*0.5;
                                j21 = (py[1]-py[-1])*0.5;
                                j31 = (pz[1]-pz[-1])*0.5;
                            }

                            if ((i1==0) || (i1==dm[1]-1))
                            {
                                 j12 = 0.0;
                                 j22 = 1.0;
                                 j32 = 0.0;
                            }
                            else
                            {
                                j12 = (px[dm[0]]-px[-dm[0]])*0.5;
                                j22 = (py[dm[0]]-py[-dm[0]])*0.5;
                                j32 = (pz[dm[0]]-pz[-dm[0]])*0.5;
                            }

                            if ((i2==0) || (i2==dm[2]-1))
                            {
                                j13 = 0.0;
                                j23 = 0.0;
                                j33 = 1.0;
                            }
                            else
                            {
                                mwSignedIndex op = dm[0]*dm[1];
                                j13 = (px[op]-px[-op])*0.5;
                                j23 = (py[op]-py[-op])*0.5;
                                j33 = (pz[op]-pz[-op])*0.5;
                            }
                        }
                    }

                    ij11 = j22*j33-j23*j32;
                    ij12 = j13*j32-j12*j33;
                    ij13 = j12*j23-j13*j22;
                    dj   = j11*ij11 + j21*ij12 + j31*ij13;
                    /* dj   = (dj*0.99+0.01); */
                    dj   = 1.0/dj;

                    ij11*= dj;
                    ij12*= dj;
                    ij13*= dj;
                    ij21 = (j23*j31-j21*j33)*dj;
                    ij22 = (j11*j33-j13*j31)*dj;
                    ij23 = (j13*j21-j11*j23)*dj;
                    ij31 = (j21*j32-j22*j31)*dj;
                    ij32 = (j12*j31-j11*j32)*dj;
                    ij33 = (j11*j22-j12*j21)*dj;

                    f1 = pf[i];
                    f2 = pf[i+my];
                    f3 = pf[i+my*2];

                    rf[0] = ij11*f1 + ij21*f2 + ij31*f3;
                    rf[1] = ij12*f1 + ij22*f2 + ij32*f3;
                    rf[2] = ij13*f1 + ij23*f2 + ij33*f3;

                    ix   = (mwSignedIndex)floor(x); dx1=x-ix; dx2=1.0-dx1;
                    iy   = (mwSignedIndex)floor(y); dy1=y-iy; dy2=1.0-dy1;
                    iz   = (mwSignedIndex)floor(z); dz1=z-iz; dz2=1.0-dz1;

                    /* Weights for trilinear interpolation */
                    w000 = dx2*dy2*dz2;
                    w100 = dx1*dy2*dz2;
                    w010 = dx2*dy1*dz2;
                    w110 = dx1*dy1*dz2;
                    w001 = dx2*dy2*dz1;
                    w101 = dx1*dy2*dz1;
                    w011 = dx2*dy1*dz1;
                    w111 = dx1*dy1*dz1;

                    ix   = bound(ix, dmo[0]);
                    iy   = bound(iy, dmo[1]);
                    iz   = bound(iz, dmo[2]);
                    ix1  = bound(ix+1, dmo[0]);
                    iy1  = bound(iy+1, dmo[1]);
                    iz1  = bound(iz+1, dmo[2]);

                    /* Neighbouring voxels used for trilinear interpolation */
                    tmpz  = dmo[1]*iz;
                    tmpy  = dmo[0]*(iy + tmpz);
                    o000  = ix +tmpy;
                    o100  = ix1+tmpy;
                    tmpy  = dmo[0]*(iy1 + tmpz);
                    o010  = ix +tmpy;
                    o110  = ix1+tmpy;
                    tmpz  = dmo[1]*iz1;
                    tmpy  = dmo[0]*(iy + tmpz);
                    o001  = ix +tmpy;
                    o101  = ix1+tmpy;
                    tmpy  = dmo[0]*(iy1 + tmpz);
                    o011  = ix +tmpy;
                    o111  = ix1+tmpy;

                    for (j=0; j<3; j++)
                    {
                        /* Increment the images themselves */
                        float *pj = po+mo*j;
                        float  f  = rf[j];
                        pj[o000] += f*w000;
                        pj[o100] += f*w100;
                        pj[o010] += f*w010;
                        pj[o110] += f*w110;
                        pj[o001] += f*w001;
                        pj[o101] += f*w101;
                        pj[o011] += f*w011;
                        pj[o111] += f*w111;
                    }
                }
            }
        }
    }
}


/*
 * t0 = Id + v0*sc
 */
void smalldef(mwSize dm[], double sc, float v0[], float t0[])
{
    mwSize j0, j1, j2;
    mwSize m = dm[0]*dm[1]*dm[2];
    float *v1 = v0+m, *v2 = v1+m;
    float *t1 = t0+m, *t2 = t1+m;

    for(j2=0; j2<dm[2]; j2++)
    {
        for(j1=0; j1<dm[1]; j1++)
        {
            for(j0=0; j0<dm[0]; j0++)
            {
                *(t0++) = (j0+1) + *(v0++)*sc;
                *(t1++) = (j1+1) + *(v1++)*sc;
                *(t2++) = (j2+1) + *(v2++)*sc;
            }
        }
    }
}

/*
 * t0 = Id + v0*sc
 * J0 = Id + I+diag(v0)*sc
 */
void smalldef_jac(mwSize dm[], double sc, float v0[], float t0[], float J0[])
{
    mwSize j0, j1, j2;
    mwSize m = dm[0]*dm[1]*dm[2];
    double sc2 = sc/2.0;
    float *v1 = v0+m, *v2 = v1+m;

    for(j2=0; j2<dm[2]; j2++)
    {
        mwSize j2m1, j2p1;
        j2m1 = bound(j2-1,dm[2]);
        j2p1 = bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSize j1m1, j1p1;
            j1m1 = bound(j1-1,dm[1]);
            j1p1 = bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSize o, om1, op1;
                o         = j0+dm[0]*(j1+dm[1]*j2);
                t0[o    ] = (j0+1) + v0[o]*sc;
                t0[o+m  ] = (j1+1) + v1[o]*sc;
                t0[o+m*2] = (j2+1) + v2[o]*sc;

                /*
                if (v0[o] > 0)
                {
                    om1 = o;
                    op1 = bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                }
                else
                {
                    om1 = bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                    op1 = o;
                }
                */

                om1 = bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                J0[o    ] = (v0[op1]-v0[om1])*sc2 + 1.0;
                J0[o+  m] = (v1[op1]-v1[om1])*sc2;
                J0[o+2*m] = (v2[op1]-v2[om1])*sc2;

                /*
                if (v1[o] > 0)
                {
                    om1 = o;
                    op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                }
                else
                {
                    om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                    op1 = o;
                }
                */

                om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                J0[o+3*m] = (v0[op1]-v0[om1])*sc2;
                J0[o+4*m] = (v1[op1]-v1[om1])*sc2 + 1.0;
                J0[o+5*m] = (v2[op1]-v2[om1])*sc2;

                /*
                if (v2[o] > 0)
                {
                    om1 = o;
                    op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                }
                {
                    om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                    op1 = o;
                }
                */

                om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                J0[o+6*m] = (v0[op1]-v0[om1])*sc2;
                J0[o+7*m] = (v1[op1]-v1[om1])*sc2;
                J0[o+8*m] = (v2[op1]-v2[om1])*sc2 + 1.0;
            }
        }
    }
}

/*
 * t0 = Id + v0*sc
 * J0 = Id + expm(D v0)
 */
void smalldef_jac1(mwSize dm[], double sc, float v0[], float t0[], float J0[])
{
    mwSize j0, j1, j2;
    mwSize m = dm[0]*dm[1]*dm[2];
    double sc2 = sc/2.0;
    float  *v1 = v0+m, *v2 = v1+m;
    float  A[9], E[9];

    for(j2=0; j2<dm[2]; j2++)
    {
        mwSize j2m1, j2p1;
        j2m1 = bound(j2-1,dm[2]);
        j2p1 = bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSize j1m1, j1p1;
            j1m1 = bound(j1-1,dm[1]);
            j1p1 = bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSize o, om1, op1;
                o         = j0+dm[0]*(j1+dm[1]*j2);
                t0[o    ] = (j0+1) + v0[o]*sc;
                t0[o+m  ] = (j1+1) + v1[o]*sc;
                t0[o+m*2] = (j2+1) + v2[o]*sc;

                om1  = bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1  = bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                A[0] = (v0[op1]-v0[om1])*sc2;
                A[1] = (v1[op1]-v1[om1])*sc2;
                A[2] = (v2[op1]-v2[om1])*sc2;

                om1  = j0+dm[0]*(j1m1+dm[1]*j2);
                op1  = j0+dm[0]*(j1p1+dm[1]*j2);
                A[3] = (v0[op1]-v0[om1])*sc2;
                A[4] = (v1[op1]-v1[om1])*sc2;
                A[5] = (v2[op1]-v2[om1])*sc2;

                om1  = j0+dm[0]*(j1+dm[1]*j2m1);
                op1  = j0+dm[0]*(j1+dm[1]*j2p1);
                A[6] = (v0[op1]-v0[om1])*sc2;
                A[7] = (v1[op1]-v1[om1])*sc2;
                A[8] = (v2[op1]-v2[om1])*sc2;

                expm22(A, E);

                J0[o    ] = E[0];
                J0[o+  m] = E[1];
                J0[o+2*m] = E[2];
                J0[o+3*m] = E[3];
                J0[o+4*m] = E[4];
                J0[o+5*m] = E[5];
                J0[o+6*m] = E[6];
                J0[o+7*m] = E[7];
                J0[o+8*m] = E[8];

            }
        }
    }
}

void divergence(mwSize dm[], float v0[], float dv[])
{
    mwSize j0, j1, j2;
    mwSize m = dm[0]*dm[1]*dm[2];
    float *v1 = v0+m, *v2 = v1+m;
    double div;

    for(j2=0; j2<dm[2]; j2++)
    {
        mwSize j2m1, j2p1;
        j2m1 = bound(j2-1,dm[2]);
        j2p1 = bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSize j1m1, j1p1;
            j1m1 = bound(j1-1,dm[1]);
            j1p1 = bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSize om1, op1;

                om1 = bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                div = v0[op1]-v0[om1];

                om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                div+= v1[op1]-v1[om1];

                om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                div+= v2[op1]-v2[om1];

                dv[j0+dm[0]*(j1+dm[1]*j2)] = 0.5*div;
            }
        }
    }
}

/* Minimum and maximum of the divergence */
void minmax_div(mwSize dm[], float v0[], double mnmx[])
{
    mwSize  j0, j1, j2;
    mwSize  m = dm[0]*dm[1]*dm[2];
    float  *v1 = v0+m, *v2 = v1+m;
    double div, maxdiv = -1e32, mindiv = 1e32;

    for(j2=0; j2<dm[2]; j2++)
    {
        mwSize j2m1, j2p1;
        j2m1 = bound(j2-1,dm[2]);
        j2p1 = bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSize j1m1, j1p1;
            j1m1 = bound(j1-1,dm[1]);
            j1p1 = bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSize om1, op1;

                om1 = bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                div = v0[op1]-v0[om1];

                om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                div+= v1[op1]-v1[om1];

                om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                div+= v2[op1]-v2[om1];

                if (div<mindiv) mindiv = div;
                if (div>maxdiv) maxdiv = div;
            }
        }
    }
    mnmx[0] = 0.5*mindiv;
    mnmx[1] = 0.5*maxdiv;
}

/*
 * Jacobian determinant field
 */
void determinant(mwSize dm[], float J0[], float d[])
{
    mwSize m = dm[0]*dm[1]*dm[2];
    mwSize j;
    double j00, j01, j02, j10, j11, j12, j20, j21, j22;

    for(j=0; j<m; j++)
    {
        j00  = J0[j    ]; j01 = J0[j+m*3]; j02 = J0[j+m*6];
        j10  = J0[j+m  ]; j11 = J0[j+m*4]; j12 = J0[j+m*7];
        j20  = J0[j+m*2]; j21 = J0[j+m*5]; j22 = J0[j+m*8];
        d[j] = j00*(j11*j22-j12*j21)+j10*(j02*j21-j01*j22)+j20*(j01*j12-j02*j11);
    }
}

/*
 * Attempt to unwrap the deformations.
 * Note: this is not always guaranteed to work,
 * but it should for most cases.
 */
void unwrap(mwSize dm[], float f[])
{
    mwSignedIndex i0, i1, i2;

   if (get_bound())
      return;

    for(i2=0; i2<dm[2]; i2++)
    {
        float *pt = f + (i2+2*dm[2])*dm[0]*dm[1];
        if (i2==0)
        {
            for(i1=0; i1<dm[1]*dm[0]; i1++)
                pt[i1] = pt[i1]-floor(pt[i1]/dm[2]+0.5)*dm[2];
        }
        else
        {
            for(i1=0; i1<dm[1]*dm[0]; i1++)
                pt[i1] = pt[i1]-floor((pt[i1]-pt[i1-dm[0]*dm[1]])/dm[2]+0.5)*dm[2];
        }
    }

    for(i1=0; i1<dm[1]; i1++)
    {
        float *pt = f + (i1+dm[2]*dm[1])*dm[0];
        if (i1==0)
        {
            for(i2=0; i2<dm[2]; i2++)
            {
                float *pt1 = pt+i2*dm[0]*dm[1];
                for(i0=0; i0<dm[0]; i0++)
                {
                    pt1[i0] = pt1[i0]-floor(pt1[i0]/dm[1]+0.5)*dm[1];
                }
            }
        }
        else
        {
            for(i2=0; i2<dm[2]; i2++)
            {
                float *pt1 = pt+i2*dm[0]*dm[1];
                for(i0=0; i0<dm[0]; i0++)
                {
                    pt1[i0] = pt1[i0]-floor((pt1[i0]-pt1[i0-dm[0]])/dm[1]+0.5)*dm[1];
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
                float *pt1 = pt + i2*dm[0]*dm[1];
                for(i1=0; i1<dm[0]*dm[1]; i1+=dm[0])
                    pt1[i1] = pt1[i1]-floor(pt1[i1]/dm[0]+0.5)*dm[0];
            }
        }
        else
        {
            for(i2=0; i2<dm[2]; i2++)
            {
                float *pt1 = pt + i2*dm[0]*dm[1];
                for(i1=0; i1<dm[0]*dm[1]; i1+=dm[0])
                    pt1[i1] = pt1[i1]-floor((pt1[i1]-pt1[i1-1])/dm[0]+0.5)*dm[0];
            }
        }
    }
}

