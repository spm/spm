/* $Id: diffeo3d.c 3802 2010-03-29 13:07:15Z john $ */
/* (c) John Ashburner (2007) */

#include <mex.h>
#include <math.h>
#include <stdio.h>
#include "optimizer3d.h"
extern double   log(double x);
extern double   exp(double x);
#define LOG(x) (((x)>0) ? log(x+0.001): -6.9078)
#define WRAP(i,m) (((i)>=0) ? (i)%(m) : ((m)+(i)%(m))%m)

/*
 * Lie Bracket
 * C = [A,B]
 */
void bracket(int dm[], float *A, float *B, float *C)
{
    float *Ax, *Ay, *Az;
    float *Bx, *By, *Bz;
    float *Cx, *Cy, *Cz;
    int i, j, k, mm = dm[0]*dm[1]*dm[2];

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
            int o1, oi1, opj1, omj1, opk1, omk1;
            o1   = dm[0]*(j+dm[1]*k);
            oi1  = dm[0]*(j+dm[1]*k);
            opj1 = dm[0]*(WRAP(j+1,dm[1])+dm[1]*k);
            omj1 = dm[0]*(WRAP(j-1,dm[1])+dm[1]*k);
            opk1 = dm[0]*(j+dm[1]*WRAP(k+1,dm[2]));
            omk1 = dm[0]*(j+dm[1]*WRAP(k-1,dm[2]));

            for(i=0; i<dm[0]; i++)
            {
                int o, opi, omi, opj, omj, opk, omk;
                double j00, j01, j02,  j10, j11, j12,  j20, j21, j22;
                double tx, ty, tz,  cx1, cy1, cz1,  cx2, cy2, cz2;

                o   = i+o1;
                opi = WRAP(i+1,dm[0])+oi1;
                omi = WRAP(i-1,dm[0])+oi1;
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

                Cx[o] = cx1-cx2;
                Cy[o] = cy1-cy2;
                Cz[o] = cz1-cz2;
            }
        }
    }
}

/*
 * In place Cholesky decomposition
 */
void chol3(int m, float A[])
{
    float *p00 = A,     *p11 = A+m,   *p22 = A+m*2,
          *p01 = A+m*3, *p02 = A+m*4, *p12 = A+m*5;
    double a00, a11, a22, a01, a02, a12;
    double s;
    int i;
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
void chol3recon(int m, float A[])
{
    float *p00 = A,     *p11 = A+m,   *p22 = A+m*2,
          *p01 = A+m*3, *p02 = A+m*4, *p12 = A+m*5;
    double a00, a11, a22, a01, a02, a12;
    int i;
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

/*
 * Composition operation
 * C(Id) = B(A(Id))
 */
void composition(int dm[], float *A, float *B, float *C)
{
    float *Ax, *Ay, *Az, *Bx, *By, *Bz, *Cx, *Cy, *Cz;
    int i, m = dm[0], n = dm[1], l = dm[2], mm = m*n*l;

    Ax = A;
    Ay = A+mm;
    Az = A+mm*2;
    Bx = B;
    By = B+mm;
    Bz = B+mm*2;
    Cx = C;
    Cy = C+mm;
    Cz = C+mm*2;

    for(i=0; i<mm; i++)
    {
        int o000,o001,o010,o011,o100,o101,o110,o111;
        double x, y, z;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        int ix, iy, iz, ix1, iy1, iz1;

        x     = Ax[i]-1.0;
        y     = Ay[i]-1.0;
        z     = Az[i]-1.0;
        ix    = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy    = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz    = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix    = WRAP(ix,m);
        iy    = WRAP(iy,n);
        iz    = WRAP(iz,l);
        ix1   = WRAP(ix+1,m);
        iy1   = WRAP(iy+1,n);
        iz1   = WRAP(iz+1,l);

        o000  = ix +m*(iy +n*iz );
        o100  = ix1+m*(iy +n*iz );
        o010  = ix +m*(iy1+n*iz );
        o110  = ix1+m*(iy1+n*iz );
        o001  = ix +m*(iy +n*iz1);
        o101  = ix1+m*(iy +n*iz1);
        o011  = ix +m*(iy1+n*iz1);
        o111  = ix1+m*(iy1+n*iz1);

        k000  = Bx[o000]-1.0;
        k100  = Bx[o100]-1.0;
        k010  = Bx[o010]-1.0;
        k110  = Bx[o110]-1.0;
        k001  = Bx[o001]-1.0;
        k101  = Bx[o101]-1.0;
        k011  = Bx[o011]-1.0;
        k111  = Bx[o111]-1.0;

        k100 -= floor((k100-k000)/m+0.5)*m;
        k010 -= floor((k010-k000)/m+0.5)*m;
        k110 -= floor((k110-k000)/m+0.5)*m;
        k001 -= floor((k001-k000)/m+0.5)*m;
        k101 -= floor((k101-k000)/m+0.5)*m;
        k011 -= floor((k011-k000)/m+0.5)*m;
        k111 -= floor((k111-k000)/m+0.5)*m;
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

        k100 -= floor((k100-k000)/l+0.5)*l;
        k010 -= floor((k010-k000)/l+0.5)*l;
        k110 -= floor((k110-k000)/l+0.5)*l;
        k001 -= floor((k001-k000)/l+0.5)*l;
        k101 -= floor((k101-k000)/l+0.5)*l;
        k011 -= floor((k011-k000)/l+0.5)*l;
        k111 -= floor((k111-k000)/l+0.5)*l;
        Cz[i] = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
              + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 + 1.0;

    }
}

/*
 * Composition operation, along with Jacobian matrices
 * C(Id)  =  B(A(Id))
 * JC(Id) = JB(A(Id))*JA(Id) ?
 */
void composition_jacobian(int dm[],
                     float *A, float *JA, float *B, float *JB,
                     float *C, float *JC)
{
    float *Ax, *Ay, *Az, *JA00, *JA01, *JA02,  *JA10, *JA11, *JA12,  *JA20, *JA21, *JA22;
    float *Bx, *By, *Bz, *JB00, *JB01, *JB02,  *JB10, *JB11, *JB12,  *JB20, *JB21, *JB22,  jb[3][3];
    float *Cx, *Cy, *Cz, *JC00, *JC01, *JC02,  *JC10, *JC11, *JC12,  *JC20, *JC21, *JC22;
    int i, mm = dm[0]*dm[1]*dm[2];

    Ax   =  A;
    Ay   =  A+mm;
    Az   =  A+mm*2;
    JA00 = JA+mm*0; JA01 = JA+mm*1; JA02 = JA+mm*2;
    JA10 = JA+mm*3; JA11 = JA+mm*4; JA12 = JA+mm*5;
    JA20 = JA+mm*6; JA21 = JA+mm*7; JA22 = JA+mm*8;

    Bx   =  B;
    By   =  B+mm;
    Bz   =  B+mm*2;
    JB00 = JB+mm*0; JB01 = JB+mm*1; JB02 = JB+mm*2;
    JB10 = JB+mm*3; JB11 = JB+mm*4; JB12 = JB+mm*5;
    JB20 = JB+mm*6; JB21 = JB+mm*7; JB22 = JB+mm*8;

    Cx   =  C;
    Cy   =  C+mm;
    Cz   =  C+mm*2;
    JC00 = JC+mm*0; JC01 = JC+mm*1; JC02 = JC+mm*2;
    JC10 = JC+mm*3; JC11 = JC+mm*4; JC12 = JC+mm*5;
    JC20 = JC+mm*6; JC21 = JC+mm*7; JC22 = JC+mm*8;

    for(i=0; i<mm; i++)
    {
        double x, y, z;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        double ja0, ja1, ja2;
        int ix, iy, iz, ix1, iy1, iz1;
        int o000, o100, o010, o110, o001, o101, o011, o111;
        int tmpz, tmpy, n;
        float *ptr;

        x    = Ax[i]-1.0;
        y    = Ay[i]-1.0;
        z    = Az[i]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = WRAP(ix  ,dm[0]);
        iy   = WRAP(iy  ,dm[1]);
        iz   = WRAP(iz  ,dm[2]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        iz1  = WRAP(iz+1,dm[2]);

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
}

/*
 * Composition operation, along with Jacobian determinants
 * C(Id)  =  B(A(Id))
 * JC(Id) = JB(A(Id))*JA(Id)
 */
void composition_jacdet(int dm[],
                     float *A, float *JA, float *B, float *JB,
                     float *C, float *JC)
{
    float *Ax, *Ay, *Az;
    float *Bx, *By, *Bz, jb;
    float *Cx, *Cy, *Cz;
    int i, mm = dm[0]*dm[1]*dm[2];

    Ax   =  A;
    Ay   =  A+mm;
    Az   =  A+mm*2;
    Bx   =  B;
    By   =  B+mm;
    Bz   =  B+mm*2;
    Cx   =  C;
    Cy   =  C+mm;
    Cz   =  C+mm*2;

    for(i=0; i<mm; i++)
    {
        double x, y, z;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        int ix, iy, iz, ix1, iy1, iz1;
        int o000, o100, o010, o110, o001, o101, o011, o111;
        int tmpz, tmpy, n;

        x    = Ax[i]-1.0;
        y    = Ay[i]-1.0;
        z    = Az[i]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = WRAP(ix  ,dm[0]);
        iy   = WRAP(iy  ,dm[1]);
        iz   = WRAP(iz  ,dm[2]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        iz1  = WRAP(iz+1,dm[2]);

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

        jb    = ((JB[o000]*dx2 + JB[o100]*dx1)*dy2 + (JB[o010]*dx2 + JB[o110]*dx1)*dy1)*dz2
              + ((JB[o001]*dx2 + JB[o101]*dx1)*dy2 + (JB[o011]*dx2 + JB[o111]*dx1)*dy1)*dz1;

        JC[i] = jb * JA[i];
    }
}

/*
 * Sample a point
 * s = f(x,y,z)
 */
double samp(int dm[], float f[], double x, double y, double z)
{
    int ix, iy, iz, ix1, iy1, iz1;
    int o000, o100, o010, o110, o001, o101, o011, o111;
    int tmpz, tmpy;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
    iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
    iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
    ix   = WRAP(ix  ,dm[0]);
    iy   = WRAP(iy  ,dm[1]);
    iz   = WRAP(iz  ,dm[2]);
    ix1  = WRAP(ix+1,dm[0]);
    iy1  = WRAP(iy+1,dm[1]);
    iz1  = WRAP(iz+1,dm[2]);

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
void sampn(int dm[], float f[], int n, int mm, double x, double y, double z, double v[])
{
    int ix, iy, iz, ix1, iy1, iz1, j;
    int o000, o100, o010, o110, o001, o101, o011, o111;
    int tmpz, tmpy;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
    iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
    iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
    ix   = WRAP(ix  ,dm[0]);
    iy   = WRAP(iy  ,dm[1]);
    iz   = WRAP(iz  ,dm[2]);
    ix1  = WRAP(ix+1,dm[0]);
    iy1  = WRAP(iy+1,dm[1]);
    iz1  = WRAP(iz+1,dm[2]);

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
void push(int dm[], int m, int n, float def[], float pf[], float po[], float so[])
{
    int   ix, iy, iz, ix1, iy1, iz1, i, j, mm;
    int   tmpz, tmpy;
    float *px, *py, *pz;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    px = def;
    py = def+m;
    pz = def+m*2;
    mm = dm[0]*dm[1]*dm[2];

    for(i=0; i<m; i++)
    {
        double x, y, z;

        if (mxIsFinite(pf[i]))
        {
            x    = px[i]-1.0; /* Subtract 1 because of MATLAB indexing */
            y    = py[i]-1.0;
            z    = pz[i]-1.0;

            /* Check range and avoid inserting values outside the FOV. */
            if (x>=1 && x<dm[0]-1 && y>=1 && y<dm[1]-1 && z>=1 && z<dm[2]-1)
            {
                /* A faster function fo voxels that are safely inside the FOV */
                int   o000, o100, o010, o110, o001, o101, o011, o111;
                float w000, w100, w010, w110, w001, w101, w011, w111;
                ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
                iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
                iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;

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
            else if ((x>=0.0) && (x<dm[0]) && (y>=0.0) && (y<dm[1]) && (z>=0.0) && (z<dm[2]))
            {
                /* A slower function for voxels at the edge of the field of view */
                int   o[8], nn=0, k;
                float w[8];

                ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
                iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
                iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
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
void pushc(int dm[], int m, int n, float def[], float pf[], float po[], float so[])
{
    int   ix, iy, iz, ix1, iy1, iz1, i, j, mm;
    int   tmpz, tmpy;
    float *px, *py, *pz;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    px = def;
    py = def+m;
    pz = def+m*2;
    mm = dm[0]*dm[1]*dm[2];

    for(i=0; i<m; i++)
    {
        double x, y, z;

        if (mxIsFinite(pf[i]))
        {
            int   o000, o100, o010, o110, o001, o101, o011, o111;
            float w000, w100, w010, w110, w001, w101, w011, w111;

            x    = px[i]-1.0; /* Subtract 1 because of MATLAB indexing */
            y    = py[i]-1.0;
            z    = pz[i]-1.0;

            ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
            iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
            iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;

            /* Weights for trilinear interpolation */
            w000 = dx2*dy2*dz2;
            w100 = dx1*dy2*dz2;
            w010 = dx2*dy1*dz2;
            w110 = dx1*dy1*dz2;
            w001 = dx2*dy2*dz1;
            w101 = dx1*dy2*dz1;
            w011 = dx2*dy1*dz1;
            w111 = dx1*dy1*dz1;

            ix   = WRAP(ix, dm[0]);
            iy   = WRAP(iy, dm[1]);
            iz   = WRAP(iz, dm[2]);
            ix1  = WRAP(ix+1, dm[0]);
            iy1  = WRAP(iy+1, dm[1]);
            iz1  = WRAP(iz+1, dm[2]);

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
void pushc_grads(int dm[], int m, float def[], float J[], float pf[], float po[])
{
    int   ix, iy, iz, ix1, iy1, iz1, i, j, mm;
    int   tmpz, tmpy;
    float *px, *py, *pz;
    float *pj11, *pj12, *pj13, *pj21, *pj22, *pj23, *pj31, *pj32, *pj33;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    px = def;
    py = def+m;
    pz = def+m*2;
    mm = dm[0]*dm[1]*dm[2];

    pj11 = J;
    pj21 = J+m;
    pj31 = J+m*2;
    pj12 = J+m*3;
    pj22 = J+m*4;
    pj32 = J+m*5;
    pj13 = J+m*6;
    pj23 = J+m*7;
    pj33 = J+m*8;

    for(i=0; i<m; i++)
    {
        double x, y, z;

        if (mxIsFinite(pf[i]))
        {
            int   o000, o100, o010, o110, o001, o101, o011, o111;
            float w000, w100, w010, w110, w001, w101, w011, w111;
            float j11, j12, j13, j21, j22, j23, j31, j32, j33;
            float ij11, ij12, ij13, ij21, ij22, ij23, ij31, ij32, ij33, dj;
            float rf[3], f1, f2, f3;

            x    = px[i]-1.0; /* Subtract 1 because of MATLAB indexing */
            y    = py[i]-1.0;
            z    = pz[i]-1.0;

            j11  = pj11[i]; j12  = pj12[i]; j13  = pj13[i];
            j21  = pj21[i]; j22  = pj22[i]; j23  = pj23[i];
            j31  = pj31[i]; j32  = pj32[i]; j33  = pj33[i];

            ij11 = j22*j33-j23*j32;
            ij12 = j13*j32-j12*j33;
            ij13 = j12*j23-j13*j22;
            dj   = j11*ij11 + j21*ij12 + j31*ij13;
            dj   = (dj*0.99+0.01);
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
            f2 = pf[i+m];
            f3 = pf[i+m*2];

            rf[0] = ij11*f1 + ij21*f2 + ij31*f3;
            rf[1] = ij12*f1 + ij22*f2 + ij32*f3;
            rf[2] = ij13*f1 + ij23*f2 + ij33*f3;

            ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
            iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
            iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;

            /* Weights for trilinear interpolation */
            w000 = dx2*dy2*dz2;
            w100 = dx1*dy2*dz2;
            w010 = dx2*dy1*dz2;
            w110 = dx1*dy1*dz2;
            w001 = dx2*dy2*dz1;
            w101 = dx1*dy2*dz1;
            w011 = dx2*dy1*dz1;
            w111 = dx1*dy1*dz1;

            ix   = WRAP(ix, dm[0]);
            iy   = WRAP(iy, dm[1]);
            iz   = WRAP(iz, dm[2]);
            ix1  = WRAP(ix+1, dm[0]);
            iy1  = WRAP(iy+1, dm[1]);
            iz1  = WRAP(iz+1, dm[2]);

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

            for (j=0; j<3; j++)
            {
                /* Increment the images themselves */
                float *pj = po+mm*j;
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


static int pow2(int k)
{
    int j0, td = 1;
    for(j0=0; j0<k; j0++)
        td = td*2;
    return(td);
}

/*
 * t0 = Id + v0*sc
 */
void smalldef(int dm[], double sc, float v0[], float t0[])
{
    int j0, j1, j2;
    int m = dm[0]*dm[1]*dm[2];
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
static void smalldef_jac(int dm[], double sc, float v0[], float t0[], float J0[])
{
    int j0, j1, j2;
    int m = dm[0]*dm[1]*dm[2];
    double sc2 = sc/2.0;
    float *v1 = v0+m, *v2 = v1+m;

    for(j2=0; j2<dm[2]; j2++)
    {
        int j2m1, j2p1;
        j2m1 = WRAP(j2-1,dm[2]);
        j2p1 = WRAP(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            int j1m1, j1p1;
            j1m1 = WRAP(j1-1,dm[1]);
            j1p1 = WRAP(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                int o, om1, op1;
                o         = j0+dm[0]*(j1+dm[1]*j2);
                t0[o    ] = (j0+1) + v0[o]*sc;
                t0[o+m  ] = (j1+1) + v1[o]*sc;
                t0[o+m*2] = (j2+1) + v2[o]*sc;

                /*
                if (v0[o] > 0)
                {
                    om1 = o;
                    op1 = WRAP(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                }
                else
                {
                    om1 = WRAP(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                    op1 = o;
                }
                */

                om1 = WRAP(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = WRAP(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
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

/* syms a[0] a[1] a[2] a[3] a[4] a[5] a[6] a[7] a[8]
   syms b[0] b[1] b[2] b[3] b[4] b[5] b[6] b[7] b[8]
   A = [a[0] a[1] a[2]; a[3] a[4] a[5]; a[6] a[7] a[8]].';
   B = [b[0] b[1] b[2]; b[3] b[4] b[5]; b[6] b[7] b[8]].';
   C = A*B;
   C = C(:)

   C = A\B;
   C = C(:)
*/
float *sub33(float *a, float *b, float *c)
{
    int i;
    for(i=0; i<9; i++)
        c[i] = a[i] - b[i];
    return(c);
}

float *add33(float *a, float *b, float *c)
{
    int i;
    for(i=0; i<9; i++)
        c[i] = a[i] + b[i];
    return(c);
}

float *scale33(float *a, float s, float *b)
{
    int i;
    for(i=0; i<9; i++)
        b[i] = a[i]*s;
    return(b);
}

float *eye33(float *a)
{
    a[0] = 1.0;
    a[1] = 0.0;
    a[2] = 0.0;
    a[3] = 0.0;
    a[4] = 1.0;
    a[5] = 0.0;
    a[6] = 0.0;
    a[7] = 0.0;
    a[8] = 1.0;
    return(a);
}

float *mul33(float *a, float *b, float *c)
{
    c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
    c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
    c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
    c[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
    c[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
    c[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
    c[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
    c[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
    c[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
    return(c);
}

float *div33(float *a, float *b, float *c)
{
    float d = a[0]*(a[4]*a[8] - a[5]*a[7]) + a[1]*(a[5]*a[6] - a[3]*a[8]) + a[2]*(a[3]*a[7] - a[4]*a[6]);
    c[0] =   (a[3]*(a[7]*b[2] - a[8]*b[1]) + a[4]*(a[8]*b[0] - a[6]*b[2]) + a[5]*(a[6]*b[1] - a[7]*b[0]))/d;
    c[1] =  -(a[0]*(a[7]*b[2] - a[8]*b[1]) + a[1]*(a[8]*b[0] - a[6]*b[2]) + a[2]*(a[6]*b[1] - a[7]*b[0]))/d;
    c[2] =   (a[0]*(a[4]*b[2] - a[5]*b[1]) + a[1]*(a[5]*b[0] - a[3]*b[2]) + a[2]*(a[3]*b[1] - a[4]*b[0]))/d;
    c[3] =   (a[3]*(a[7]*b[5] - a[8]*b[4]) + a[4]*(a[8]*b[3] - a[6]*b[5]) + a[5]*(a[6]*b[4] - a[7]*b[3]))/d;
    c[4] =  -(a[0]*(a[7]*b[5] - a[8]*b[4]) + a[1]*(a[8]*b[3] - a[6]*b[5]) + a[2]*(a[6]*b[4] - a[7]*b[3]))/d;
    c[5] =   (a[0]*(a[4]*b[5] - a[5]*b[4]) + a[1]*(a[5]*b[3] - a[3]*b[5]) + a[2]*(a[3]*b[4] - a[4]*b[3]))/d;
    c[6] =   (a[3]*(a[7]*b[8] - a[8]*b[7]) + a[4]*(a[8]*b[6] - a[6]*b[8]) + a[5]*(a[6]*b[7] - a[7]*b[6]))/d;
    c[7] =  -(a[0]*(a[7]*b[8] - a[8]*b[7]) + a[1]*(a[8]*b[6] - a[6]*b[8]) + a[2]*(a[6]*b[7] - a[7]*b[6]))/d;
    c[8] =   (a[0]*(a[4]*b[8] - a[5]*b[7]) + a[1]*(a[5]*b[6] - a[3]*b[8]) + a[2]*(a[3]*b[7] - a[4]*b[6]))/d;
    return(c);
}

void pade33(float *a, float *l)
{
    float u[9], v[9], num[9], den[9], a0[9], a2[9], a3[9];
    
    eye33(a0);
    mul33(a,a,a2);
    mul33(a2,a,a3);
    scale33(a0,120.0,a0);
    scale33(a2, 12.0,a2);
    add33(a0,a2,u);
    scale33(a,60.0,v);
    add33(v,a3,v);
    
    add33(u,v,num);
    sub33(u,v,den);
    div33(den,num,l);
}

void pade22(float *a, float *l)
{
    float u[9], v[9], num[9], den[9], a0[9], a2[9];
    
    eye33(a0);
    mul33(a,a,a2);
    scale33(a0,12.0,a0);
    add33(a0,a2,u);
    scale33(a,6.0,v);
    
    add33(u,v,num);
    sub33(u,v,den);
    div33(den,num,l);
}

float norm1(float *a)
{
    float r, rm;
    rm = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);
    r  = fabs(a[3]) + fabs(a[4]) + fabs(a[5]);
    if (r>rm) rm = r;
    r  = fabs(a[6]) + fabs(a[7]) + fabs(a[8]);
    if (r>rm) rm = r;
    return(rm);
}

float *assign33(float *a, float *b)
{
    int i;
    for(i=0; i<9; i++)
        b[i] = a[i];
    return(b);
}

void expm33(float *a, float *l)
{
    /* See expm.m in MATLAB or http://mathworld.wolfram.com/PadeApproximant.html */
    int K;
    K = (int)ceil(log((double)(norm1(a)*2.3481))*1.44269504088896);
    if (K>0)
    {
        float b[9];
        float s = 1.0/pow2(K);
        int i;
        scale33(a,s,b);
        pade33(b, l);
        for(i=0; i<K; i++)
        {
            assign33(l,b);
            mul33(b,b,l);
        }
    }
    else
        pade33(a, l);
}

void expm22(float *a, float *l)
{
    /* See expm.m in MATLAB or http://mathworld.wolfram.com/PadeApproximant.html */
    int K;
    K = (int)ceil(log((double)(norm1(a)*12.356))*1.44269504088896);
    if (K>0)
    {
        float b[9];
        float s = 1.0/pow2(K);
        int i;
        scale33(a,s,b);
        pade22(b, l);
        for(i=0; i<K; i++)
        {
            assign33(l,b);
            mul33(b,b,l);
        }
    }
    else
        pade22(a, l);
}

/*
 * t0 = Id + v0*sc
 * J0 = Id + expm(D v0)
 */
void smalldef_jac1(int dm[], double sc, float v0[], float t0[], float J0[])
{
    int j0, j1, j2;
    int m = dm[0]*dm[1]*dm[2];
    double sc2 = sc/2.0;
    float *v1 = v0+m, *v2 = v1+m;
    float A[9], E[9];

    for(j2=0; j2<dm[2]; j2++)
    {
        int j2m1, j2p1;
        j2m1 = WRAP(j2-1,dm[2]);
        j2p1 = WRAP(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            int j1m1, j1p1;
            j1m1 = WRAP(j1-1,dm[1]);
            j1p1 = WRAP(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                int o, om1, op1;
                o         = j0+dm[0]*(j1+dm[1]*j2);
                t0[o    ] = (j0+1) + v0[o]*sc;
                t0[o+m  ] = (j1+1) + v1[o]*sc;
                t0[o+m*2] = (j2+1) + v2[o]*sc;

                om1  = WRAP(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1  = WRAP(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
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

/* Minimum and maximum of the divergence */
void minmax_div(int dm[], float v0[], double mnmx[])
{
    int j0, j1, j2;
    int m = dm[0]*dm[1]*dm[2];
    float *v1 = v0+m, *v2 = v1+m;
    float div, maxdiv = -1e32, mindiv = 1e32;

    for(j2=0; j2<dm[2]; j2++)
    {
        int j2m1, j2p1;
        j2m1 = WRAP(j2-1,dm[2]);
        j2p1 = WRAP(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            int j1m1, j1p1;
            j1m1 = WRAP(j1-1,dm[1]);
            j1p1 = WRAP(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                int om1, op1;

                om1 = WRAP(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = WRAP(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
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
    mnmx[0] = (double)mindiv;
    mnmx[1] = (double)maxdiv;
}

/*
 * t0 = Id + v0*sc
 * J0 = Id + |I+diag(v0)*sc|
 */
static void smalldef_jacdet(int dm[], double sc, float v0[], float t0[], float J0[])
{
    int j0, j1, j2;
    int m = dm[0]*dm[1]*dm[2];
    double sc2 = sc/2.0;
    float *v1 = v0+m, *v2 = v1+m;
    
    for(j2=0; j2<dm[2]; j2++)
    {
        int j2m1, j2p1;
        j2m1 = WRAP(j2-1,dm[2]);
        j2p1 = WRAP(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            int j1m1, j1p1;
            j1m1 = WRAP(j1-1,dm[1]);
            j1p1 = WRAP(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                int o, om1, op1;
                double j00,j10,j20, j01,j11,j21, j02,j12,j22;
                
                o         = j0+dm[0]*(j1+dm[1]*j2);
                t0[o    ] = (j0+1) + v0[o]*sc;
                t0[o+m  ] = (j1+1) + v1[o]*sc;
                t0[o+m*2] = (j2+1) + v2[o]*sc;

                om1 = WRAP(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = WRAP(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                
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
 * Jacobian determinant field
 */
void determinant(int dm[], float J0[], float d[])
{
    int m = dm[0]*dm[1]*dm[2];
    int j;
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
 * J0 := J0*inv(I+diag(v0)*sc)
 */
static void jac_div_smalldef(int dm[], double sc, float v0[], float J0[])
{
    int j0, j1, j2;
    int m = dm[0]*dm[1]*dm[2];
    double sc2 = sc/2.0;
    float *v1 = v0+m, *v2 = v1+m;

    for(j2=0; j2<dm[2]; j2++)
    {
        int j2m1, j2p1;
        j2m1 = WRAP(j2-1,dm[2]);
        j2p1 = WRAP(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            int j1m1, j1p1;
            j1m1 = WRAP(j1-1,dm[1]);
            j1p1 = WRAP(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                int o, om1, op1;
                double j00,j01,j02, j10,j11,j12, j20,j21,j22;
                double t00,t01,t02, t10,t11,t12, t20,t21,t22;
                double idt;

                om1 = WRAP(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = WRAP(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
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
void expdef(int dm[], int k, double sc, float v[], float t0[], float t1[], float J0[], float J1[])
{
    float *optr;
    int m = dm[0]*dm[1]*dm[2];
    int j;

    optr = t0;

    if(J0!=(float *)0)
    {
        smalldef_jac(dm, sc/pow2(k), v, t0, J0);
        for(j=0; j<k; j++)
        {
            float *tmpp;
            composition_jacobian(dm, t0, J0, t0, J0, t1, J1);
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
            composition(dm, t0, t0, t1);
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
void expdefdet(int dm[], int k, double sc, float v[], float t0[], float t1[], float J0[], float J1[])
{
    float *optr;
    int m = dm[0]*dm[1]*dm[2];
    int j;

    optr = t0;

    if(J0!=(float *)0)
    {
        smalldef_jacdet(dm, sc/pow2(k), v, t0, J0);
        for(j=0; j<k; j++)
        {
            float *tmpp;
            composition_jacdet(dm, t0, J0, t0, J0, t1, J1);
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
            composition(dm, t0, t0, t1);
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


static double smalldef_objfun_mn(int dm[], float f[], float g[], float v[], float jd[], double sc, float b[], float A[])
{
    int j, j0, j1, j2, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;

    j = 0;
    for(j2=0; j2<dm[2]; j2++)
        for(j1=0; j1<dm[1]; j1++)
            for(j0=0; j0<dm[0]; j0++, j++)
    {
        double x, y, z;
        int    ix, iy, iz, ix1, iy1, iz1, k;
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
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = WRAP(ix,dm[0]);
        iy   = WRAP(iy,dm[1]);
        iz   = WRAP(iz,dm[2]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        iz1  = WRAP(iz+1,dm[2]);
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
                int km = k*m;
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
                int k1;
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

static double smalldef_objfun2(int dm[], float f[], float g[], float v[], float jd[], double sc, float b[], float A[])
{
    int j, j0, j1, j2, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;

    j = 0;
    for(j2=0; j2<dm[2]; j2++)
        for(j1=0; j1<dm[1]; j1++)
            for(j0=0; j0<dm[0]; j0++, j++)
    {
        double x, y, z;
        int   ix, iy, iz, ix1, iy1, iz1, k;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        double d, dx, dy, dz, sd, sdx, sdy, sdz, ss;
        
        x    = j0 + sc*v[j    ];
        y    = j1 + sc*v[j+m  ];
        z    = j2 + sc*v[j+m*2];
        
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = WRAP(ix,dm[0]);
        iy   = WRAP(iy,dm[1]);
        iz   = WRAP(iz,dm[2]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        iz1  = WRAP(iz+1,dm[2]);

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
            int km = k*m;
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

static double smalldef_objfun(int dm[], float f[], float g[], float v[], float jd[], double sc, float b[], float A[])
{
    int j,j0,j1,j2, m = dm[0]*dm[1]*dm[2];
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
        int   ix, iy, iz, ix1, iy1, iz1;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx1, dx2, dy1, dy2, dz1, dz2;
        double d, dx, dy, dz, dt = 1.0;
        
        x    = j0 + sc*v[j    ];
        y    = j1 + sc*v[j+m  ];
        z    = j2 + sc*v[j+m*2];

        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = WRAP(ix,dm[0]);
        iy   = WRAP(iy,dm[1]);
        iz   = WRAP(iz,dm[2]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        iz1  = WRAP(iz+1,dm[2]);
        
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

static double initialise_objfun_mn(int dm[], float f[], float g[], float t0[], float J0[], float jd[], float b[], float A[])
{
    int j, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;
 
    for(j=0; j<m; j++)
    {
        double x, y, z;
        int    ix, iy, iz, ix1, iy1, iz1, k;
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
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = WRAP(ix,dm[0]);
        iy   = WRAP(iy,dm[1]);
        iz   = WRAP(iz,dm[2]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        iz1  = WRAP(iz+1,dm[2]);
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
                int km = k*m;
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
                int k1;
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

static double initialise_objfun2(int dm[], float f[], float g[], float t0[], float J0[], float jd[], float b[], float A[])
{
    int j, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0;
    
    for(j=0; j<m; j++)
    {
        double x, y, z;
        int   ix, iy, iz, ix1, iy1, iz1, k;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx0, dx1, dx2, dy0, dy1, dy2, dz0, dz1, dz2;
        double d, dx, dy, dz;
        double sd, sdx, sdy, sdz, ss;

        x    = t0[j    ]-1.0;
        y    = t0[j+m  ]-1.0;
        z    = t0[j+m*2]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = WRAP(ix,dm[0]);
        iy   = WRAP(iy,dm[1]);
        iz   = WRAP(iz,dm[2]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        iz1  = WRAP(iz+1,dm[2]);

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
            int km = k*m;
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

static double initialise_objfun(int dm[], float f[], float g[], float t0[], float J0[], float jd[], float b[], float A[])
{
    int j, m = dm[0]*dm[1]*dm[2];
    double ssl = 0.0, dt = 1.0;

    if (dm[3]>1)
    {
        return(initialise_objfun2(dm, f, g, t0, J0, jd, b, A));
    }

    for(j=0; j<m; j++)
    {
        double x, y, z;
        int   ix, iy, iz, ix1, iy1, iz1;
        double k000, k100, k010, k110, k001, k101, k011, k111;
        double dx0, dx1, dx2, dy0, dy1, dy2, dz0, dz1, dz2;
        double d, dx, dy, dz;

        x    = t0[j    ]-1.0;
        y    = t0[j+m  ]-1.0;
        z    = t0[j+m*2]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        iz   = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;
        ix   = WRAP(ix,dm[0]);
        iy   = WRAP(iy,dm[1]);
        iz   = WRAP(iz,dm[2]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        iz1  = WRAP(iz+1,dm[2]);

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


static void squaring(int dm[], int k, int save_transf, float b[], float A[], float t0[], float t1[], float J0[], float J1[])
{
    int i, j, m = dm[0]*dm[1]*dm[2];
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
            sampn(dm, b, 3, m, x, y, z, bs);
            b0 = bs[0];
            b1 = bs[1];
            b2 = bs[2];

            buf1[j    ] = dt*(b0*j00+b1*j10+b2*j20);
            buf1[j+m  ] = dt*(b0*j01+b1*j11+b2*j21);
            buf1[j+m*2] = dt*(b0*j02+b1*j12+b2*j22);

            /* J'*A*J */
            sampn(dm, A, 6, m, x, y, z, as);
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
            composition_jacobian(dm, t0, J0, t0, J0, t1, J1);
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

/*
 * Attempt to unwrap the deformations.
 * Note: this is not always guarunteed to work,
 * but it should for most cases.
 */
void unwrap(int dm[], float f[])
{
    int i0, i1, i2;

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

int iteration_scratchsize(int dm[], int code, int k)
{
    int m1, m2;
    int m = dm[0]*dm[1]*dm[2];
    if (k>0)
    {
        m1 = 30*m;
        if (code==1) m1 += 9*m;
        m2 = 9*m+fmg3_scratchsize(dm);
        if (m1>m2)
            return(m1);
        else
            return(m2);
    }
    else
    {
        m1 = 9*m;
        if (code==1) m1 += 6*m;
        m2 = 9*m + fmg3_scratchsize(dm);
        if (m1>m2)
            return(m1);
        else
            return(m2);
    }
}

void iteration(int dm[], int k, float v[], float g[], float f[], float jd[],
               int rtype, double param0[], double lmreg0, int cycles, int its, int code,
               float ov[], double ll[], float *buf)
{
    float *sbuf;
    float *b, *A;
    double ssl, ssp, sc;
    static double param[6] = {1.0,1.0,1.0,1.0,0.0,0.0};
    int m = dm[0]*dm[1]*dm[2];
    int j;

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

    if (rtype==0)
        LtLf_le(dm, v, param, sbuf);
    else if (rtype==1)
        LtLf_me(dm, v, param, sbuf);
    else /* if (rtype==2) */
        LtLf_be(dm, v, param, sbuf);

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

    if (lmreg0>0.0) param[5] = param[5] + lmreg0;

    fmg3(dm, A, b, rtype, param, cycles, its, sbuf, sbuf+3*m); 
    for(j=0; j<m*3; j++) ov[j] = v[j] - sbuf[j];
}

