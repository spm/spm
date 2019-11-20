/* $Id: shoot_diffeo3d.c 7699 2019-11-20 16:00:03Z yael $ */
/* (c) John Ashburner (2011) */

#include <math.h>
#include <stdio.h>
#include "spm_mex.h"
#include "shoot_optim3d.h"
#include "shoot_expm3.h"
#include "shoot_boundary.h"

extern double   log(double x);
extern double   exp(double x);
#define LOG(x) (((x)>0) ? log(x+0.001): -6.9078)


static void grad1(mwSize d0, mwSize d1, mwSize d2, mwSize N, int boundary_code, float F[], float D[])
{
    mwSize k, n;
    float *fp = F, *dp = D, *Dend;

    if (d1>1)
    {
        for(n=0; n<N; n++)
        {
            for(k=0;k<d2; k++)
            {
                /* element 0 */
                switch(boundary_code)
                {
                    case BOUND_CIRCULANT:
                        for(Dend = dp+d0; dp<Dend; *(dp++) = (*(fp+d0) - *(fp+d0*(d1-1)))/2.0f, fp++);
                        break;
                    case BOUND_NEUMANN:
                        for(Dend = dp+d0; dp<Dend; *(dp++) = (*(fp+d0) - *fp)/2.0f, fp++);
                        break;
                    case BOUND_DIRICHLET:
                        for(Dend = dp+d0; dp<Dend; *(dp++) = (*(fp+d0))/2.0f, fp++);
                        break;
                }

                /* elements 1..d1-2 */
                for(Dend = dp+d0*(d1-2); dp<Dend; *(dp++) = (*(fp+d0) - *(fp-d0))/2.0f, fp++);

                /* element d1-1 */
                switch(boundary_code)
                {
                    case BOUND_CIRCULANT:
                        for(Dend = dp+d0; dp<Dend; *(dp++) = (*(fp-d0*(d1-1)) - *(fp-d0))/2.0f, fp++);
                        break;
                    case BOUND_NEUMANN:
                        for(Dend = dp+d0; dp<Dend; *(dp++) = (*fp - *(fp-d0))/2.0f, fp++);
                        break;
                    case BOUND_DIRICHLET:
                        for(Dend = dp+d0; dp<Dend; *(dp++) = -(*(fp-d0))/2.0f, fp++);
                        break;
                }
            }
        }
    }
    else
        for(Dend = D+N*d2*d1*d0; dp<Dend; *(dp++) = 0.0);
}

void grad(mwSize dm[], mwSize N, float F[], float D[])
{
    int boundary_code = get_bound();
    mwSize M = N*dm[0]*dm[1]*dm[2];
    grad1(1,           dm[0], dm[1]*dm[2], N, boundary_code, F, D+M*0); /* x-gradient */
    grad1(dm[0],       dm[1],       dm[2], N, boundary_code, F, D+M*1); /* y-gradient */
    grad1(dm[0]*dm[1], dm[2],           1, N, boundary_code, F, D+M*2); /* z-gradient */
}



/*
 * Lie Bracket
 * C = [A,B]
 *
 * See https://en.wikipedia.org/wiki/Lie_bracket_of_vector_fields
 *
 * This code is not used anywhere, but retained for reference in
 * case it is needed in future. Note that effects at the boundaries
 * may need fixing.
*/
void bracket(mwSize dm[], float *A, float *B, float *C)
{
    float *A1, *A2, *A3;
    float *B1, *B2, *B3;
    float *C1, *C2, *C3;
    mwSize i, j, k, mm = dm[0]*dm[1]*dm[2];

    /* Components of first vector field */
    A1 = A;
    A2 = A + mm;
    A3 = A + mm*2;

    /* Components of 2nd vector field */
    B1 = B;
    B2 = B + mm;
    B3 = B + mm*2;

    /* Components of returned vector field */
    C1 = C;
    C2 = C + mm;
    C3 = C + mm*2;

    for(k=0; k<dm[2]; k++) /* loop over z dimension */
    {
        for(j=0; j<dm[1]; j++) /* loop over y dimension */
        {
            mwSize o1, opj1, omj1, opk1, omk1;

            /* Compute some offsets for later use */
            o1   = dm[0]*(j+dm[1]*k);
            opj1 = dm[0]*(bound(j+1,dm[1])+dm[1]*k);
            omj1 = dm[0]*(bound(j-1,dm[1])+dm[1]*k);
            opk1 = dm[0]*(j+dm[1]*bound(k+1,dm[2]));
            omk1 = dm[0]*(j+dm[1]*bound(k-1,dm[2]));

            for(i=0; i<dm[0]; i++) /* loop over x dimension */
            {
                mwSize o, opi, omi, opj, omj, opk, omk;
                float j00, j01, j02,  j10, j11, j12,  j20, j21, j22;
                float tx, ty, tz,  cx1, cy1, cz1,  cx2, cy2, cz2;

                /* Neighbour indices for computing Jacobian matrices */
                o   = i+o1;
                opi = (mwSize)bound(i+1,dm[0])+o1;
                omi = (mwSize)bound(i-1,dm[0])+o1;
                opj = i+opj1;
                omj = i+omj1;
                opk = i+opk1;
                omk = i+omk1;

                /* x, y and z components of A deformation */
                tx = A1[o];
                ty = A2[o];
                tz = A3[o];

                /* Elements of Jacobian matrix of B deformation */
                j00 = (B1[opi]-B1[omi])/2.0f;
                j01 = (B2[opi]-B2[omi])/2.0f;
                j02 = (B3[opi]-B3[omi])/2.0f;

                j10 = (B1[opj]-B1[omj])/2.0f;
                j11 = (B2[opj]-B2[omj])/2.0f;
                j12 = (B3[opj]-B3[omj])/2.0f;

                j20 = (B1[opk]-B1[omk])/2.0f;
                j21 = (B2[opk]-B2[omk])/2.0f;
                j22 = (B3[opk]-B3[omk])/2.0f;

                /* Multiply transforms by Jacobian */
                cx1 = tx*j00+ty*j10+tz*j20;
                cy1 = tx*j01+ty*j11+tz*j21;
                cz1 = tx*j02+ty*j12+tz*j22;

                /* x, y and z components of B deformation */
                tx = B1[o];
                ty = B2[o];
                tz = B3[o];

                /* Elements of Jacobian matrix of A deformation */
                j00 = (A1[opi]-A1[omi])/2.0f;
                j01 = (A2[opi]-A2[omi])/2.0f;
                j02 = (A3[opi]-A3[omi])/2.0f;

                j10 = (A1[opj]-A1[omj])/2.0f;
                j11 = (A2[opj]-A2[omj])/2.0f;
                j12 = (A3[opj]-A3[omj])/2.0f;

                j20 = (A1[opk]-A1[omk])/2.0f;
                j21 = (A2[opk]-A2[omk])/2.0f;
                j22 = (A3[opk]-A3[omk])/2.0f;

                /* Multiply transforms by Jacobian */
                cx2 = tx*j00+ty*j10+tz*j20;
                cy2 = tx*j01+ty*j11+tz*j21;
                cz2 = tx*j02+ty*j12+tz*j22;

                /* Compute Lie bracket */
                C1[o] = cx2-cx1;
                C2[o] = cy2-cy1;
                C3[o] = cz2-cz1;
            }
        }
    }
}

/*
 * Composition operations, possibly along with Jacobian matrices.
 * Done using B(A) = U(A) + A, where B = Id+U.
 */
static void composition_stuff(mwSize dm[], mwSize mm,
                              float *B, /*@null@*/ float *JB, float *A, /*@null@*/ float *JA,
                              float *C, /*@null@*/ float *JC, int flag)
{
    float *A1, *A2, *A3, *JA00, *JA01, *JA02,  *JA10, *JA11, *JA12,  *JA20, *JA21, *JA22;
    float *B1, *B2, *B3, *JB00, *JB01, *JB02,  *JB10, *JB11, *JB12,  *JB20, *JB21, *JB22;
    float *C1, *C2, *C3, *JC00, *JC01, *JC02,  *JC10, *JC11, *JC12,  *JC20, *JC21, *JC22;
    float  NaN = (float)mxGetNaN();
    mwSize i, mmb = dm[0]*dm[1]*dm[2];

    A1   =  A;
    A2   =  A+mm;
    A3   =  A+mm*2;
    B1   =  B;
    B2   =  B+mmb;
    B3   =  B+mmb*2;
    C1   =  C;
    C2   =  C+mm;
    C3   =  C+mm*2;

    /* Only relevant if (JC!=0 && flag==0)
     * Note that tensors are transposed, which was not deliberate.
     * Maybe refactor this to add some clarity. */
    JA00 = JA+mm*0; JA01 = JA+mm*1; JA02 = JA+mm*2;
    JA10 = JA+mm*3; JA11 = JA+mm*4; JA12 = JA+mm*5;
    JA20 = JA+mm*6; JA21 = JA+mm*7; JA22 = JA+mm*8;

    JB00 = JB+mmb*0; JB01 = JB+mmb*1; JB02 = JB+mmb*2;
    JB10 = JB+mmb*3; JB11 = JB+mmb*4; JB12 = JB+mmb*5;
    JB20 = JB+mmb*6; JB21 = JB+mmb*7; JB22 = JB+mmb*8;

    JC00 = JC+mm*0; JC01 = JC+mm*1; JC02 = JC+mm*2;
    JC10 = JC+mm*3; JC11 = JC+mm*4; JC12 = JC+mm*5;
    JC20 = JC+mm*6; JC21 = JC+mm*7; JC22 = JC+mm*8;

    for(i=0; i<mm; i++)
    {
        float x, y, z;
        float k000, k100, k010, k110, k001, k101, k011, k111;
        float dx1, dx2, dy1, dy2, dz1, dz2;
        mwSignedIndex ix, iy, iz, ix0, iy0, iz0, ix1, iy1, iz1;
        mwSize o000, o100, o010, o110, o001, o101, o011, o111;
        mwSize tmpz, tmpy;

        x    = A1[i]-1.0f;
        y    = A2[i]-1.0f;
        z    = A3[i]-1.0f;

        if (mxIsFinite(x) && mxIsFinite(y) && mxIsFinite(z))
        {
            ix   = (mwSignedIndex)floor((double)x); dx1=x-(float)ix; dx2=1.0f-dx1;
            iy   = (mwSignedIndex)floor((double)y); dy1=y-(float)iy; dy2=1.0f-dy1;
            iz   = (mwSignedIndex)floor((double)z); dz1=z-(float)iz; dz2=1.0f-dz1;
            ix0  = bound(ix  ,dm[0]);
            iy0  = bound(iy  ,dm[1]);
            iz0  = bound(iz  ,dm[2]);
            ix1  = bound(ix+1,dm[0]);
            iy1  = bound(iy+1,dm[1]);
            iz1  = bound(iz+1,dm[2]);

            tmpz  = dm[1]*iz0;
            tmpy  = dm[0]*(iy0 + tmpz);
            o000  = (mwSize)(ix0 +tmpy);
            o100  = (mwSize)(ix1+tmpy);
            tmpy  = dm[0]*(iy1 + tmpz);
            o010  = (mwSize)(ix0 +tmpy);
            o110  = (mwSize)(ix1+tmpy);
            tmpz  = dm[1]*iz1;
            tmpy  = dm[0]*(iy0 + tmpz);
            o001  = (mwSize)(ix0 +tmpy);
            o101  = (mwSize)(ix1+tmpy);
            tmpy  = dm[0]*(iy1 + tmpz);
            o011  = (mwSize)(ix0 +tmpy);
            o111  = (mwSize)(ix1+tmpy);

            /* Sample from x-displacement field - subtract x coordinates of sampled points */
            k000  = B1[o000] - (float)(ix0+1);
            k100  = B1[o100] - (float)(ix1+1);
            k010  = B1[o010] - (float)(ix0+1);
            k110  = B1[o110] - (float)(ix1+1);
            k001  = B1[o001] - (float)(ix0+1);
            k101  = B1[o101] - (float)(ix1+1);
            k011  = B1[o011] - (float)(ix0+1);
            k111  = B1[o111] - (float)(ix1+1);
            C1[i] = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                  + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 + x + 1.0f;

            /* sample from y-displacement field */
            k000  = B2[o000] - (float)(iy0+1);
            k100  = B2[o100] - (float)(iy0+1);
            k010  = B2[o010] - (float)(iy1+1);
            k110  = B2[o110] - (float)(iy1+1);
            k001  = B2[o001] - (float)(iy0+1);
            k101  = B2[o101] - (float)(iy0+1);
            k011  = B2[o011] - (float)(iy1+1);
            k111  = B2[o111] - (float)(iy1+1);
            C2[i] = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                  + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 + y + 1.0f;

            /* sample from z-displacement field */
            k000  = B3[o000] - (float)(iz0+1);
            k100  = B3[o100] - (float)(iz0+1);
            k010  = B3[o010] - (float)(iz0+1);
            k110  = B3[o110] - (float)(iz0+1);
            k001  = B3[o001] - (float)(iz1+1);
            k101  = B3[o101] - (float)(iz1+1);
            k011  = B3[o011] - (float)(iz1+1);
            k111  = B3[o111] - (float)(iz1+1);
            C3[i] = ((k000*dx2 + k100*dx1)*dy2 + (k010*dx2 + k110*dx1)*dy1)*dz2
                  + ((k001*dx2 + k101*dx1)*dy2 + (k011*dx2 + k111*dx1)*dy1)*dz1 + z + 1.0f;

            if (JC!=0 && JB!=0 && JA!=0)
            {
                if (flag==0)
                {
                    float *ptr;
                    float ja0, ja1, ja2;
                    float jb[3][3];

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
                    float jb;
                    jb    = ((JB[o000]*dx2 + JB[o100]*dx1)*dy2 + (JB[o010]*dx2 + JB[o110]*dx1)*dy1)*dz2
                          + ((JB[o001]*dx2 + JB[o101]*dx1)*dy2 + (JB[o011]*dx2 + JB[o111]*dx1)*dy1)*dz1;
                    JC[i] = jb * JA[i];
                }
            }
        }
        else
        {
            C1[i] = NaN;
            C2[i] = NaN;
            C3[i] = NaN;
            if (JC!=0 && JB!=0 && JA!=0)
            {
                if (flag==0)
                {
                    JC00[i] = NaN;
                    JC01[i] = NaN;
                    JC02[i] = NaN;
                    JC10[i] = NaN;
                    JC11[i] = NaN;
                    JC12[i] = NaN;
                    JC20[i] = NaN;
                    JC21[i] = NaN;
                    JC22[i] = NaN;
                }
                else
                {
                    JC[i] = NaN;
                }
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
    composition_stuff(dm, mm, B, 0, A, 0, C, 0, 0);
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

static void def2jac_wrap(mwSize dm[], float *Psi, float *Jpsi, mwSignedIndex s, int code)
{
    mwSize i, j, k, k0, k2, mm;

    if (s!=-1)
    {
        k0 = (mwSize)(s-1);
        k2 = (mwSize)s;
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
            y1 = Psi + dm[0]*(j+dm[1]*k);
            y2 = Psi + dm[0]*(j+dm[1]*(k+dm[2]));
            y3 = Psi + dm[0]*(j+dm[1]*(k+dm[2]*2));
            dp = Jpsi + dm[0]*(j+dm[1]*k);

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im, ip;
                float j11, j12, j13, j21, j22, j23, j31, j32, j33;
                im = bound(i-1,dm[0])-i;
                ip = bound(i+1,dm[0])-i;

                j11 = (y1[i+ip]-y1[i+im])/(float)dm[0]; j11 = 0.5f*(j11 - (float)floor(j11+0.5))*(float)dm[0];
                j21 = (y2[i+ip]-y2[i+im])/(float)dm[0]; j21 = 0.5f*(j21 - (float)floor(j21+0.5))*(float)dm[0];
                j31 = (y3[i+ip]-y3[i+im])/(float)dm[0]; j31 = 0.5f*(j31 - (float)floor(j31+0.5))*(float)dm[0];

                j12 = (y1[i+jp]-y1[i+jm])/(float)dm[1]; j12 = 0.5f*(j12 - (float)floor(j12+0.5))*(float)dm[1];
                j22 = (y2[i+jp]-y2[i+jm])/(float)dm[1]; j22 = 0.5f*(j22 - (float)floor(j22+0.5))*(float)dm[1];
                j32 = (y3[i+jp]-y3[i+jm])/(float)dm[1]; j32 = 0.5f*(j32 - (float)floor(j32+0.5))*(float)dm[1];

                if (dm[2]>1)
                {
                    j13 = (y1[i+kp]-y1[i+km])/(float)dm[2]; j13 = 0.5f*(j13 - (float)floor(j13+0.5))*(float)dm[2];
                    j23 = (y2[i+kp]-y2[i+km])/(float)dm[2]; j23 = 0.5f*(j23 - (float)floor(j23+0.5))*(float)dm[2];
                    j33 = (y3[i+kp]-y3[i+km])/(float)dm[2]; j33 = 0.5f*(j33 - (float)floor(j33+0.5))*(float)dm[2];
                }
                else
                {
                    j13 = 0.0f;
                    j23 = 0.0f;
                    j33 = 1.0f;
                }

                if (code==0)
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

static void def2jac_neuman(mwSize dm[], float *Psi, float *Jpsi, mwSignedIndex s, int code)
{
    mwSize mm, k, k0, k2;

    if (s!=-1)
    {
        k0 = (mwSize)(s-1);
        k2 = (mwSize)s;
    }
    else
    {
        k0 = 0;
        k2 = dm[2];
    }
    mm = dm[0]*dm[1]*(k2-k0);

    for(k=k0; k<k2; k++)
    {
        mwSize j;
        for(j=0; j<dm[1]; j++)
        {
            mwSize i;
            float *y1, *y2, *y3, *dp;
            y1 = Psi + dm[0]*(j+dm[1]*k);
            y2 = Psi + dm[0]*(j+dm[1]*(k+dm[2]));
            y3 = Psi + dm[0]*(j+dm[1]*(k+dm[2]*2));
            dp = Jpsi + dm[0]*(j+dm[1]*k);

            for(i=0; i<dm[0]; i++)
            {
                float j11, j12, j13, j21, j22, j23, j31, j32, j33;

                if ((i==0) || (i==dm[0]-1))
                {
                    j11 = 1.0f;
                    j21 = 0.0f;
                    j31 = 0.0f;
                }
                else
                {
                    j11 = (y1[i+1]-y1[i-1])*0.5f;
                    j21 = (y2[i+1]-y2[i-1])*0.5f;
                    j31 = (y3[i+1]-y3[i-1])*0.5f;
                }

                if ((j==0) || (j==dm[1]-1))
                {
                    j12 = 0.0f;
                    j22 = 1.0f;
                    j32 = 0.0f;
                }
                else
                {
                    mwSignedIndex op = (mwSignedIndex)i+dm[0], om = (mwSignedIndex)i-dm[0];
                    j12 = (y1[op]-y1[om])*0.5f;
                    j22 = (y2[op]-y2[om])*0.5f;
                    j32 = (y3[op]-y3[om])*0.5f;
                }

                if ((k==0) || (k==(dm[2])-1))
                {
                    j13 = 0.0f;
                    j23 = 0.0f;
                    j33 = 1.0f;
                }
                else
                {
                    mwSignedIndex op = (mwSignedIndex)i+dm[0]*dm[1], om = (mwSignedIndex)i-dm[0]*dm[1];
                    j13 = (y1[op]-y1[om])*0.5f;
                    j23 = (y2[op]-y2[om])*0.5f;
                    j33 = (y3[op]-y3[om])*0.5f;
                }

                if (code==0)
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

void def2det(mwSize dm[], float *Psi, float *Jpsi, mwSignedIndex s)
{
    if (get_bound()!=0)
        def2jac_neuman(dm, Psi, Jpsi, s, 0);
    else
        def2jac_wrap(dm, Psi, Jpsi, s, 0);
}

void def2jac(mwSize dm[], float *Psi, float *Jpsi, mwSignedIndex s)
{
    if (get_bound()!=0)
        def2jac_neuman(dm, Psi, Jpsi, s, 1);
    else
        def2jac_wrap(dm, Psi, Jpsi, s, 1);
}

#define TINY 5e-2f

void pushpull(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F0[], /*@null@@out@*/float S0[], float F1[], unsigned int code)
{
    mwSize i, m0;
    float  *px, *py, *pz;
    float  NaN = (float)mxGetNaN();

    px = Psi;
    py = Psi+m1;
    pz = Psi+m1*2;
    m0 = dm0[0]*dm0[1]*dm0[2];

#   pragma omp parallel for
    for (i=0; i<m1; i++)
    {
        float  x, y, z;
        mwSize j;

        x    = px[i]-1.0f; /* Subtract 1 because of MATLAB indexing */
        y    = py[i]-1.0f;
        z    = pz[i]-1.0f;

	if (((code & 2)==2 && mxIsFinite(x) && mxIsFinite(y) && mxIsFinite(z))
	                   || ((x>=-TINY) && (x<=(float)(dm0[0])-1.0f+TINY)
	                   &&  (y>=-TINY) && (y<=(float)(dm0[1])-1.0f+TINY)
			   &&  (z>=-TINY) && (z<=(float)(dm0[2])-1.0f+TINY)))
        {
            mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
            mwSize o000, o100, o010, o110, o001, o101, o011, o111;
            mwSize tmpz, tmpy;
            float  dx1, dx2, dy1, dy2, dz1, dz2;

            ix   = (mwSignedIndex)floor((double)x); dx1=x-(float)ix; dx2=1.0f-dx1;
            iy   = (mwSignedIndex)floor((double)y); dy1=y-(float)iy; dy2=1.0f-dy1;
            iz   = (mwSignedIndex)floor((double)z); dz1=z-(float)iz; dz2=1.0f-dz1;
            ix   = bound(ix  ,dm0[0]);
            iy   = bound(iy  ,dm0[1]);
            iz   = bound(iz  ,dm0[2]);
            ix1  = bound(ix+1,dm0[0]);
            iy1  = bound(iy+1,dm0[1]);
            iz1  = bound(iz+1,dm0[2]);

            tmpz  = dm0[1] * iz;
            tmpy  = dm0[0] *(iy+tmpz);
            o000  = (mwSize)(ix +tmpy);
            o100  = (mwSize)(ix1+tmpy);
            tmpy  = dm0[0] *(iy1+tmpz);
            o010  = (mwSize)(ix +tmpy);
            o110  = (mwSize)(ix1+tmpy);
            tmpz  = dm0[1] * iz1;
            tmpy  = dm0[0] *(iy +tmpz);
            o001  = (mwSize)(ix +tmpy);
            o101  = (mwSize)(ix1+tmpy);
            tmpy  = dm0[0] *(iy1+tmpz);
            o011  = (mwSize)(ix +tmpy);
            o111  = (mwSize)(ix1+tmpy);

            if ((code&1)==1)
            {
                float *pf0 = F0;
                for(j=0; j<n; j++, pf0 += m0)
                {
                    F1[i+j*m1] = ((pf0[o000]*dx2 + pf0[o100]*dx1)*dy2 + (pf0[o010]*dx2 + pf0[o110]*dx1)*dy1)*dz2
                               + ((pf0[o001]*dx2 + pf0[o101]*dx1)*dy2 + (pf0[o011]*dx2 + pf0[o111]*dx1)*dy1)*dz1;
                }
            }
            else
            {
                float  w000, w100, w010, w110, w001, w101, w011, w111;
                float *pf0=F0;

                /* Weights for trilinear interpolation */
                w000 = dx2*dy2*dz2;
                w100 = dx1*dy2*dz2;
                w010 = dx2*dy1*dz2;
                w110 = dx1*dy1*dz2;
                w001 = dx2*dy2*dz1;
                w101 = dx1*dy2*dz1;
                w011 = dx2*dy1*dz1;
                w111 = dx1*dy1*dz1;

                for (j=0; j<n; j++, pf0+=m0)
                {
                    /* Increment the images themselves */
                    float  f1 = F1[i+j*m1];
                    if (mxIsFinite((double)f1))
                    {
#                       pragma omp atomic
                        pf0[o000] += f1*w000;
#                       pragma omp atomic
                        pf0[o100] += f1*w100;
#                       pragma omp atomic
                        pf0[o010] += f1*w010;
#                       pragma omp atomic
                        pf0[o110] += f1*w110;
#                       pragma omp atomic
                        pf0[o001] += f1*w001;
#                       pragma omp atomic
                        pf0[o101] += f1*w101;
#                       pragma omp atomic
                        pf0[o011] += f1*w011;
#                       pragma omp atomic
                        pf0[o111] += f1*w111;
                    }
                }

                if (S0!=0)
                {
                    /* Increment an image containing the number of voxels added
                       - based on finite values in the 1st image */
                    if (mxIsFinite((double)F1[i]))
                    {
#                       pragma omp atomic
                        S0[o000] += w000;
#                       pragma omp atomic
                        S0[o100] += w100;
#                       pragma omp atomic
                        S0[o010] += w010;
#                       pragma omp atomic
                        S0[o110] += w110;
#                       pragma omp atomic
                        S0[o001] += w001;
#                       pragma omp atomic
                        S0[o101] += w101;
#                       pragma omp atomic
                        S0[o011] += w011;
#                       pragma omp atomic
                        S0[o111] += w111;
                    }
                }
            }
        }
        else if ((code&1)==1)
        {
            for(j=0; j<n; j++)
            {
                F1[i+j*m1] = NaN;
            }
        }
    }
}

void pullc(mwSize dm0[], mwSize m1, mwSize n, float Psi[],  float F0[], /*@out@*/ float F1[])
{
    pushpull(dm0, m1, n, Psi, F0, 0, F1, 3);
}

void  pull(mwSize dm0[], mwSize m1, mwSize n, float Psi[],  float F0[], /*@out@*/ float F1[])
{
    pushpull(dm0, m1, n, Psi, F0, 0, F1, 1);
}


/* Rather than sample from an image according to a deformation,
 * it is also possible to push voxels from one image into
 * another according to the inverse of the deformation.
 * Note that the result is a noisy version of a Jacobian "modulated"
 * image.
 */
void pushc(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[], /*@out@*/ float F0[], /*@null@@out@*/ float S0[])
{
    pushpull(dm0, m1, n, Psi, F0, S0, F1, 2);
}

void  push(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[], /*@out@*/ float F0[], /*@null@@out@*/ float S0[])
{
    pushpull(dm0, m1, n, Psi, F0, S0, F1, 0);
}


/* This should perhaps be tidied up and generalised in line with the pushpull function */
void pushc_grads(mwSize dmo[], mwSize dm[], float Psi[], float Jpsi[], float pf[], float po[])
{
    mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
    mwSize i2, i, mo, my;
    float *px, *py, *pz;
    float *pj11, *pj12, *pj13, *pj21, *pj22, *pj23, *pj31, *pj32, *pj33;

    my = dm[0]*dm[1]*dm[2];
    px = Psi;
    py = Psi+my;
    pz = Psi+my*2;
    mo = dmo[0]*dmo[1]*dmo[2];

    /* Only used if Jpsi!=(float *)0 */
    pj11 = Jpsi;
    pj21 = Jpsi+my;
    pj31 = Jpsi+my*2;
    pj12 = Jpsi+my*3;
    pj22 = Jpsi+my*4;
    pj32 = Jpsi+my*5;
    pj13 = Jpsi+my*6;
    pj23 = Jpsi+my*7;
    pj33 = Jpsi+my*8;

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
                    float x, y, z;
                    float rf[3];
                    float dx1, dx2, dy1, dy2, dz1, dz2;

                    x    = *px-1.0f; /* Subtract 1 because of MATLAB indexing */
                    y    = *py-1.0f;
                    z    = *pz-1.0f;

                    if (Jpsi!=(float *)0)
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

                            j11 = (px[i0p]-px[i0m])/(float)dm[0]; j11 = 0.5f*(j11 - (float)floor(j11+0.5))*(float)dm[0];
                            j21 = (py[i0p]-py[i0m])/(float)dm[0]; j21 = 0.5f*(j21 - (float)floor(j21+0.5))*(float)dm[0];
                            j31 = (pz[i0p]-pz[i0m])/(float)dm[0]; j31 = 0.5f*(j31 - (float)floor(j31+0.5))*(float)dm[0];

                            j12 = (px[i1p]-px[i1m])/(float)dm[1]; j12 = 0.5f*(j12 - (float)floor(j12+0.5))*(float)dm[1];
                            j22 = (py[i1p]-py[i1m])/(float)dm[1]; j22 = 0.5f*(j22 - (float)floor(j22+0.5))*(float)dm[1];
                            j32 = (pz[i1p]-pz[i1m])/(float)dm[1]; j32 = 0.5f*(j32 - (float)floor(j32+0.5))*(float)dm[1];

                            if (dm[2]>1)
                            {
                                j13 = (px[i2p]-px[i2m])/(float)dm[2]; j13 = 0.5f*(j13 - (float)floor(j13+0.5))*(float)dm[2];
                                j23 = (py[i2p]-py[i2m])/(float)dm[2]; j23 = 0.5f*(j23 - (float)floor(j23+0.5))*(float)dm[2];
                                j33 = (pz[i2p]-pz[i2m])/(float)dm[2]; j33 = 0.5f*(j33 - (float)floor(j33+0.5))*(float)dm[2];
                            }
                            else
                            {
                                j13 = 0.0f;
                                j23 = 0.0f;
                                j33 = 1.0f;
                            }
                        }
                        else
                        {
                            /* Neumann boundary - only discontinuities at edges */
                            if ((i0==0) || (i0==dm[0]-1))
                            {
                                j11 = 1.0f;
                                j21 = 0.0f;
                                j31 = 0.0f;
                            }
                            else
                            {
                                j11 = (px[1]-px[-1])*0.5f;
                                j21 = (py[1]-py[-1])*0.5f;
                                j31 = (pz[1]-pz[-1])*0.5f;
                            }

                            if ((i1==0) || (i1==dm[1]-1))
                            {
                                 j12 = 0.0f;
                                 j22 = 1.0f;
                                 j32 = 0.0f;
                            }
                            else
                            {
                                j12 = (px[dm[0]]-px[-dm[0]])*0.5f;
                                j22 = (py[dm[0]]-py[-dm[0]])*0.5f;
                                j32 = (pz[dm[0]]-pz[-dm[0]])*0.5f;
                            }

                            if ((i2==0) || (i2==dm[2]-1))
                            {
                                j13 = 0.0f;
                                j23 = 0.0f;
                                j33 = 1.0f;
                            }
                            else
                            {
                                mwSignedIndex op = (mwSignedIndex)(dm[0]*dm[1]);
                                j13 = (px[op]-px[-op])*0.5f;
                                j23 = (py[op]-py[-op])*0.5f;
                                j33 = (pz[op]-pz[-op])*0.5f;
                            }
                        }
                    }

                    ij11 = j22*j33-j23*j32;
                    ij12 = j13*j32-j12*j33;
                    ij13 = j12*j23-j13*j22;
                    dj   = j11*ij11 + j21*ij12 + j31*ij13;
                    /* dj   = (dj*0.99+0.01); */
                    dj   = 1.0f/dj;

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

                    ix   = (mwSignedIndex)floor((double)x); dx1=x-(float)ix; dx2=1.0f-dx1;
                    iy   = (mwSignedIndex)floor((double)y); dy1=y-(float)iy; dy2=1.0f-dy1;
                    iz   = (mwSignedIndex)floor((double)z); dz1=z-(float)iz; dz2=1.0f-dz1;

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
                    o000  = (mwSize)(ix +tmpy);
                    o100  = (mwSize)(ix1+tmpy);
                    tmpy  = dmo[0]*(iy1 + tmpz);
                    o010  = (mwSize)(ix +tmpy);
                    o110  = (mwSize)(ix1+tmpy);
                    tmpz  = dmo[1]*iz1;
                    tmpy  = dmo[0]*(iy + tmpz);
                    o001  = (mwSize)(ix +tmpy);
                    o101  = (mwSize)(ix1+tmpy);
                    tmpy  = dmo[0]*(iy1 + tmpz);
                    o011  = (mwSize)(ix +tmpy);
                    o111  = (mwSize)(ix1+tmpy);

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
 * Psi0 = Id + V0 sc
 */
void smalldef(mwSize dm[], float sc, float V[], float Psi[])
{
    mwSize j0, j1, j2;
    mwSize m = dm[0]*dm[1]*dm[2];
    float *V0 = V, *V1 = V+m, *V2 = V+2*m;
    float *Psi0 = Psi, *Psi1 = Psi+m, *Psi2 = Psi+2*m;

    for(j2=0; j2<dm[2]; j2++)
    {
        for(j1=0; j1<dm[1]; j1++)
        {
            for(j0=0; j0<dm[0]; j0++)
            {
                *(Psi0++) = (float)(j0+1) + *(V0++)*sc;
                *(Psi1++) = (float)(j1+1) + *(V1++)*sc;
                *(Psi2++) = (float)(j2+1) + *(V2++)*sc;
            }
        }
    }
}

/*
 * Psi0  = Id + V0 sc
 * Jpsi0 = Id + I+D V0 sc
 */
void smalldef_jac(mwSize dm[], float sc, float V[], float Psi[], float Jpsi[])
{
    mwSize j0, j1, j2;
    mwSize m = dm[0]*dm[1]*dm[2];
    float sc2 = (float)(sc/2.0);
    float *V0 = V, *V1 = V+m, *V2 = V+2*m;

    for(j2=0; j2<dm[2]; j2++)
    {
        mwSize j2m1, j2p1;
        j2m1 = (mwSize)bound(j2-1,dm[2]);
        j2p1 = (mwSize)bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSize j1m1, j1p1;
            j1m1 = (mwSize)bound(j1-1,dm[1]);
            j1p1 = (mwSize)bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSize o, om1, op1;
                o         = j0+dm[0]*(j1+dm[1]*j2);
                Psi[o    ] = (float)(j0+1) + V0[o]*sc;
                Psi[o+m  ] = (float)(j1+1) + V1[o]*sc;
                Psi[o+m*2] = (float)(j2+1) + V2[o]*sc;

                om1 = (mwSize)bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = (mwSize)bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                Jpsi[o    ] = (V0[op1]-V0[om1])*sc2 + 1.0f;
                Jpsi[o+  m] = (V1[op1]-V1[om1])*sc2;
                Jpsi[o+2*m] = (V2[op1]-V2[om1])*sc2;

                om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                Jpsi[o+3*m] = (V0[op1]-V0[om1])*sc2;
                Jpsi[o+4*m] = (V1[op1]-V1[om1])*sc2 + 1.0f;
                Jpsi[o+5*m] = (V2[op1]-V2[om1])*sc2;

                om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                Jpsi[o+6*m] = (V0[op1]-V0[om1])*sc2;
                Jpsi[o+7*m] = (V1[op1]-V1[om1])*sc2;
                Jpsi[o+8*m] = (V2[op1]-V2[om1])*sc2 + 1.0f;
            }
        }
    }
}

/*
 * Psi0  = Id + V0 sc
 * Jpsi0 = Id + expm(D V0 sc)
 */
void smalldef_jac1(mwSize dm[], float sc, float V[], float Psi[], float Jpsi[])
{
    mwSize j0, j1, j2;
    mwSize m = dm[0]*dm[1]*dm[2];
    float sc2 = (float)(sc/2.0);
    float *V0 = V, *V1 = V+m, *V2 = V+2*m;
    float  A[9], E[9];

    for(j2=0; j2<dm[2]; j2++)
    {
        mwSize j2m1, j2p1;
        j2m1 = (mwSize)bound(j2-1,dm[2]);
        j2p1 = (mwSize)bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSize j1m1, j1p1;
            j1m1 = (mwSize)bound(j1-1,dm[1]);
            j1p1 = (mwSize)bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSize o, om1, op1;
                o         = j0+dm[0]*(j1+dm[1]*j2);
                Psi[o    ] = (float)(j0+1) + V0[o]*sc;
                Psi[o+m  ] = (float)(j1+1) + V1[o]*sc;
                Psi[o+m*2] = (float)(j2+1) + V2[o]*sc;

                om1  = (mwSize)bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1  = (mwSize)bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                A[0] = (V0[op1]-V0[om1])*sc2;
                A[1] = (V1[op1]-V1[om1])*sc2;
                A[2] = (V2[op1]-V2[om1])*sc2;

                om1  = j0+dm[0]*(j1m1+dm[1]*j2);
                op1  = j0+dm[0]*(j1p1+dm[1]*j2);
                A[3] = (V0[op1]-V0[om1])*sc2;
                A[4] = (V1[op1]-V1[om1])*sc2;
                A[5] = (V2[op1]-V2[om1])*sc2;

                om1  = j0+dm[0]*(j1+dm[1]*j2m1);
                op1  = j0+dm[0]*(j1+dm[1]*j2p1);
                A[6] = (V0[op1]-V0[om1])*sc2;
                A[7] = (V1[op1]-V1[om1])*sc2;
                A[8] = (V2[op1]-V2[om1])*sc2;

                expm22(A, E);

                Jpsi[o    ] = E[0];
                Jpsi[o+  m] = E[1];
                Jpsi[o+2*m] = E[2];
                Jpsi[o+3*m] = E[3];
                Jpsi[o+4*m] = E[4];
                Jpsi[o+5*m] = E[5];
                Jpsi[o+6*m] = E[6];
                Jpsi[o+7*m] = E[7];
                Jpsi[o+8*m] = E[8];

            }
        }
    }
}

/*
 * Compute divergence of a vector field.
 * dv = div(V), where V is of size dm[0]*dm[1]*dm[2]*3
 * See https://en.wikipedia.org/wiki/Divergence
 */
void divergence(mwSize dm[], float V[], float dv[])
{
    mwSize j0, j1, j2;
    mwSize m = dm[0]*dm[1]*dm[2];
    float *V0 = V, *V1 = V+m, *V2 = V+2*m; /* x, y and z components of V */
    float div;

    for(j2=0; j2<dm[2]; j2++) /* loop over z dimension */
    {
        mwSize j2m1, j2p1;
        j2m1 = (mwSize)bound(j2-1,dm[2]);
        j2p1 = (mwSize)bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++) /* loop over y dimension */
        {
            mwSize j1m1, j1p1;
            j1m1 = (mwSize)bound(j1-1,dm[1]);
            j1p1 = (mwSize)bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++) /* loop over x dimension */
            {
                mwSize om1, op1;

                /* x gradient (times 2) */
                om1 = (mwSize)bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = (mwSize)bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                div = V0[op1]-V0[om1];

                /* y gradient (times 2) */
                om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                div+= V1[op1]-V1[om1];

                /* z gradient (times 2) */
                om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                div+= V2[op1]-V2[om1];

                /* Divergence - need to half the computed sum */
                dv[j0+dm[0]*(j1+dm[1]*j2)] = 0.5f*div;
            }
        }
    }
}

/* Minimum and maximum of the divergence */
void minmax_div(mwSize dm[], float V[], double mnmx[])
{
    mwSize  j0, j1, j2;
    mwSize  m = dm[0]*dm[1]*dm[2];
    float  *V0 = V, *V1 = V+m, *V2 = V+2*m;
    double div, maxdiv = -1e32, mindiv = 1e32;

    for(j2=0; j2<dm[2]; j2++)
    {
        mwSize j2m1, j2p1;
        j2m1 = (mwSize)bound(j2-1,dm[2]);
        j2p1 = (mwSize)bound(j2+1,dm[2]);

        for(j1=0; j1<dm[1]; j1++)
        {
            mwSize j1m1, j1p1;
            j1m1 = (mwSize)bound(j1-1,dm[1]);
            j1p1 = (mwSize)bound(j1+1,dm[1]);

            for(j0=0; j0<dm[0]; j0++)
            {
                mwSize om1, op1;

                om1 = (mwSize)bound(j0-1,dm[0])+dm[0]*(j1+dm[1]*j2);
                op1 = (mwSize)bound(j0+1,dm[0])+dm[0]*(j1+dm[1]*j2);
                div = V0[op1]-V0[om1];

                om1 = j0+dm[0]*(j1m1+dm[1]*j2);
                op1 = j0+dm[0]*(j1p1+dm[1]*j2);
                div+= V1[op1]-V1[om1];

                om1 = j0+dm[0]*(j1+dm[1]*j2m1);
                op1 = j0+dm[0]*(j1+dm[1]*j2p1);
                div+= V2[op1]-V2[om1];

                if (div<mindiv) mindiv = div;
                if (div>maxdiv) maxdiv = div;
            }
        }
    }
    mnmx[0] = 0.5*mindiv;
    mnmx[1] = 0.5*maxdiv;
}

/*
 * Jacobian determinant field from a field of Jacobian matrices (dm[0]*dm[1]*dm[2]*3*3).
 * d = det(J)
 * See https://en.wikipedia.org/wiki/Determinant
 */
void determinant(mwSize dm[], float Jpsi[], float d[])
{
    mwSize m = dm[0]*dm[1]*dm[2];
    mwSize j;
    double j00, j01, j02, j10, j11, j12, j20, j21, j22;

    for(j=0; j<m; j++)
    {
        j00  = Jpsi[j    ]; j01 = Jpsi[j+m*3]; j02 = Jpsi[j+m*6];
        j10  = Jpsi[j+m  ]; j11 = Jpsi[j+m*4]; j12 = Jpsi[j+m*7];
        j20  = Jpsi[j+m*2]; j21 = Jpsi[j+m*5]; j22 = Jpsi[j+m*8];
        d[j] = (float)(j00*(j11*j22-j12*j21)+j10*(j02*j21-j01*j22)+j20*(j01*j12-j02*j11));
    }
}

