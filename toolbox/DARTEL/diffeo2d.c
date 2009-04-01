/* $Id: diffeo2d.c 3032 2009-04-01 14:14:18Z guillaume $ */
/* (c) John Ashburner (2007) */

#include <math.h>
#include<mex.h>
#include "optimizer2d.h"

#define WRAP(i,m) (((i)>=0) ? (i)%(m) : ((m)+(i)%(m))%m)

void bracket(int dm[], double *A, double *B, double *C)
{
    double *Ax, *Ay;
    double *Bx, *By;
    double *Cx, *Cy;
    int i, j;

    Ax   =   A;
    Ay   =  &A[dm[0]*dm[1]];

    Bx   =   B;
    By   =  &B[dm[0]*dm[1]];

    Cx   =   C;
    Cy   =  &C[dm[0]*dm[1]];

    for(j=0; j<dm[1]; j++)
    {
        for(i=0; i<dm[0]; i++)
        {
            int o, opi, omi, opj, omj;
            double j00, j01, j10, j11;
            double tx, ty, cx1, cy1, cx2, cy2;

            o   = i+dm[0]*j;
            opi = WRAP(i+1,dm[0])+dm[0]*j;
            omi = WRAP(i-1,dm[0])+dm[0]*j;
            opj = i+dm[0]*WRAP(j+1,dm[1]);
            omj = i+dm[0]*WRAP(j-1,dm[1]);

            tx = Ax[o];
            ty = Ay[o];
            j00 = (Bx[opi]-Bx[omi])/2.0;
            j01 = (By[opi]-By[omi])/2.0;
            j10 = (Bx[opj]-Bx[omj])/2.0;
            j11 = (By[opj]-By[omj])/2.0;
            cx1 = tx*j00+ty*j10;
            cy1 = tx*j01+ty*j11;

            tx = Bx[o];
            ty = By[o];
            j00 = (Ax[opi]-Ax[omi])/2.0;
            j01 = (Ay[opi]-Ay[omi])/2.0;
            j10 = (Ax[opj]-Ax[omj])/2.0;
            j11 = (Ay[opj]-Ay[omj])/2.0;

            cx2 = tx*j00+ty*j10;
            cy2 = tx*j01+ty*j11;

            Cx[o] = cx1-cx2;
            Cy[o] = cy1-cy2;
        }
    }
}

void composition(int dm[], double *A, double *B, double *C)
{
    double *Ax, *Ay, *Bx, *By, *Cx, *Cy;
    int i, m = dm[0], n = dm[1], mm = m*n;

    Ax = A;
    Ay = &A[mm];
    Bx = B;
    By = &B[mm];
    Cx = C;
    Cy = &C[mm];

    for(i=0; i<mm; i++)
    {
        double x, y;
        double k11,k12,k21,k22;
        double dx1, dx2, dy1, dy2;
        int ix, iy, ix1, iy1;

        x     = Ax[i]-1.0;
        y     = Ay[i]-1.0;
        ix    = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy    = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        ix    = WRAP(ix,m);
        iy    = WRAP(iy,n);
        ix1   = WRAP(ix+1,m);
        iy1   = WRAP(iy+1,n);

        k22   = Bx[ix  + m*iy ] - 1.0;
        k12   = Bx[ix1 + m*iy ] - 1.0;
        k21   = Bx[ix  + m*iy1] - 1.0;
        k11   = Bx[ix1 + m*iy1] - 1.0;
        k12   = k12-floor((k12-k22)/m+0.5)*m;
        k21   = k21-floor((k21-k22)/m+0.5)*m;
        k11   = k11-floor((k11-k22)/m+0.5)*m;
        Cx[i] = (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22   = By[ix  + m*iy ] - 1.0;
        k12   = By[ix1 + m*iy ] - 1.0;
        k21   = By[ix  + m*iy1] - 1.0;
        k11   = By[ix1 + m*iy1] - 1.0;
        k12   = k12-floor((k12-k22)/n+0.5)*n;
        k21   = k21-floor((k21-k22)/n+0.5)*n;
        k11   = k11-floor((k11-k22)/n+0.5)*n;
        Cy[i] = (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;
    }
}

void composition_jacobian(int dm[],
                     double *A, double *JA, double *B, double *JB,
                     double *C, double *JC)
{
    double *Ax, *Ay, *JA00, *JA01, *JA10, *JA11;
    double *Bx, *By, *JB00, *JB01, *JB10, *JB11, jb[2][2];
    double *Cx, *Cy, *JC00, *JC01, *JC10, *JC11;
    int i, m = dm[0], n = dm[1], mm = m*n;
    Ax   =   A;
    Ay   =  &A[mm];
    JA00 = &JA[mm*0];
    JA01 = &JA[mm*1];
    JA10 = &JA[mm*2];
    JA11 = &JA[mm*3];

    Bx   =   B;
    By   =  &B[mm];
    JB00 = &JB[mm*0];
    JB01 = &JB[mm*1];
    JB10 = &JB[mm*2];
    JB11 = &JB[mm*3];

    Cx   =   C;
    Cy   =  &C[mm];
    JC00 = &JC[mm*0];
    JC01 = &JC[mm*1];
    JC10 = &JC[mm*2];
    JC11 = &JC[mm*3];

    for(i=0; i<mm; i++)
    {
        double x, y;
        double k11,k12,k21,k22;
        double dx1, dx2, dy1, dy2;
        int ix, iy, ix1, iy1, o11,o12,o21,o22;

        x    = Ax[i]-1.0;
        y    = Ay[i]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        ix   = WRAP(ix,m);
        iy   = WRAP(iy,n);
        ix1  = WRAP(ix+1,m);
        iy1  = WRAP(iy+1,n);
        o22  = ix  + m*iy;
        o12  = ix1 + m*iy;
        o21  = ix  + m*iy1;
        o11  = ix1 + m*iy1;

        k22  = Bx[o22]-1.0;
        k12  = Bx[o12]-1.0;
        k21  = Bx[o21]-1.0;
        k11  = Bx[o11]-1.0;
        k12  = k12-floor((k12-k22)/m+0.5)*m;
        k21  = k21-floor((k21-k22)/m+0.5)*m;
        k11  = k11-floor((k11-k22)/m+0.5)*m;
        Cx[i]= (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22  = By[o22]-1.0;
        k12  = By[o12]-1.0;
        k21  = By[o21]-1.0;
        k11  = By[o11]-1.0;
        k12  = k12-floor((k12-k22)/n+0.5)*n;
        k21  = k21-floor((k21-k22)/n+0.5)*n;
        k11  = k11-floor((k11-k22)/n+0.5)*n;
        Cy[i]= (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22  = JB00[o22];
        k12  = JB00[o12];
        k21  = JB00[o21];
        k11  = JB00[o11];
        jb[0][0] =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        k22  = JB01[o22];
        k12  = JB01[o12];
        k21  = JB01[o21];
        k11  = JB01[o11];
        jb[0][1] =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        k22  = JB10[o22];
        k12  = JB10[o12];
        k21  = JB10[o21];
        k11  = JB10[o11];
        jb[1][0] =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        k22  = JB11[o22];
        k12  = JB11[o12];
        k21  = JB11[o21];
        k11  = JB11[o11];
        jb[1][1] =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        JC00[i] = jb[0][0]*JA00[i] + jb[1][0]*JA01[i];
        JC01[i] = jb[0][1]*JA00[i] + jb[1][1]*JA01[i];
        JC10[i] = jb[0][0]*JA10[i] + jb[1][0]*JA11[i];
        JC11[i] = jb[0][1]*JA10[i] + jb[1][1]*JA11[i];
    }
}

void composition_detjac(int dm[],
                     double *A, double *dA, double *B, double *dB,
                     double *C, double *dC)
{
    double *Ax, *Ay;
    double *Bx, *By, jb;
    double *Cx, *Cy;
    int i, m = dm[0], n = dm[1], mm = m*n;
    Ax   =   A;
    Ay   =  &A[mm];

    Bx   =   B;
    By   =  &B[mm];

    Cx   =   C;
    Cy   =  &C[mm];

    for(i=0; i<mm; i++)
    {
        double x, y;
        double k11,k12,k21,k22;
        double dx1, dx2, dy1, dy2;
        int ix, iy, ix1, iy1, o11,o12,o21,o22;

        x    = Ax[i]-1.0;
        y    = Ay[i]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        ix   = WRAP(ix,m);
        iy   = WRAP(iy,n);
        ix1  = WRAP(ix+1,m);
        iy1  = WRAP(iy+1,n);
        o22  = ix  + m*iy;
        o12  = ix1 + m*iy;
        o21  = ix  + m*iy1;
        o11  = ix1 + m*iy1;

        k22  = Bx[o22]-1.0;
        k12  = Bx[o12]-1.0;
        k21  = Bx[o21]-1.0;
        k11  = Bx[o11]-1.0;
        k12  = k12-floor((k12-k22)/m+0.5)*m;
        k21  = k21-floor((k21-k22)/m+0.5)*m;
        k11  = k11-floor((k11-k22)/m+0.5)*m;
        Cx[i]= (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22  = By[o22]-1.0;
        k12  = By[o12]-1.0;
        k21  = By[o21]-1.0;
        k11  = By[o11]-1.0;
        k12  = k12-floor((k12-k22)/n+0.5)*n;
        k21  = k21-floor((k21-k22)/n+0.5)*n;
        k11  = k11-floor((k11-k22)/n+0.5)*n;
        Cy[i]= (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22  = dB[o22];
        k12  = dB[o12];
        k21  = dB[o21];
        k11  = dB[o11];
        jb   =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        dC[i] = jb*dA[i];
    }
}


double samp(int dm[], double f[], double x, double y)
{
    int ix, iy, ix1, iy1, o11, o12, o21, o22;
    double dx1, dx2, dy1, dy2;

    ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
    iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
    ix   = WRAP(ix,dm[0]);
    iy   = WRAP(iy,dm[1]);
    ix1  = WRAP(ix+1,dm[0]);
    iy1  = WRAP(iy+1,dm[1]);
    o22  = ix  + dm[0]*iy;
    o12  = ix1 + dm[0]*iy;
    o21  = ix  + dm[0]*iy1;
    o11  = ix1 + dm[0]*iy1;
    return((f[o22]*dx2 + f[o12]*dx1)*dy2 + (f[o21]*dx2 + f[o11]*dx1)*dy1);
}

void expdef(int dm[], int k, double v[], double t0[], double t1[], double J0[], double J1[])
{
    double *optr;
    double td;
    int m = dm[0]*dm[1];
    int i, j;

    optr = t0;

    td = 1;
    for(i=0; i<k; i++)
        td = td*2;
    td = 1.0/td;

    if(J0!=(double *)0)
    {
        for(j=0; j<dm[1]; j++)
        {
            for(i=0; i<dm[0]; i++)
            {
                int o   = i+dm[0]* j;
                t0[o  ] = (i+1) + v[o  ]*td;
                t0[o+m] = (j+1) + v[o+m]*td;
                J0[o    ] = (v[WRAP(i+1,dm[0])+dm[0]*j  ]-v[WRAP(i-1,dm[0])+dm[0]*j  ])*td/2 + 1.0;
                J0[o+  m] = (v[WRAP(i+1,dm[0])+dm[0]*j+m]-v[WRAP(i-1,dm[0])+dm[0]*j+m])*td/2;
                J0[o+2*m] = (v[i+dm[0]*WRAP(j+1,dm[1])  ]-v[i+dm[0]*WRAP(j-1,dm[1])  ])*td/2;
                J0[o+3*m] = (v[i+dm[0]*WRAP(j+1,dm[1])+m]-v[i+dm[0]*WRAP(j-1,dm[1])+m])*td/2 + 1.0;
            }
        }
        for(i=0; i<k; i++)
        {
            double *tmpp;
            composition_jacobian(dm, t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    }
    else
    {
        for(j=0; j<dm[1]; j++)
        {
            for(i=0; i<dm[0]; i++)
            {
                int o   = i+dm[0]* j;
                t0[o  ] = (i+1) + v[o  ]*td;
                t0[o+m] = (j+1) + v[o+m]*td;
            }
        }
        for(i=0; i<k; i++)
        {
            double *tmpp;
            composition(dm, t0, t0, t1);
            tmpp = t0; t0   = t1; t1   = tmpp;
        }
    }

    if (optr != t0)
    {
        for(i=0; i<2*m; i++)
            t1[i] = t0[i];

        if (J0!=(double *)0)
            for(i=0; i<4*m; i++)
                J1[i] = J0[i];
    }
}

void expdefdet(int dm[], int k, double v[], double t0[], double t1[], double J0[], double J1[])
{
    double *optr;
    double td;
    int m = dm[0]*dm[1];
    int i, j;

    optr = t0;

    td = 1;
    for(i=0; i<k; i++)
        td = td*2;
    td = 1.0/td;

    if(J0!=(double *)0)
    {
        for(j=0; j<dm[1]; j++)
        {
            for(i=0; i<dm[0]; i++)
            {
                double j00,j01,j10,j11;
                int o   = i+dm[0]* j;
                t0[o  ] = (i+1) + v[o  ]*td;
                t0[o+m] = (j+1) + v[o+m]*td;
                j00     = (v[WRAP(i+1,dm[0])+dm[0]*j  ]-v[WRAP(i-1,dm[0])+dm[0]*j  ])*td/2 + 1.0;
                j10     = (v[WRAP(i+1,dm[0])+dm[0]*j+m]-v[WRAP(i-1,dm[0])+dm[0]*j+m])*td/2;
                j01     = (v[i+dm[0]*WRAP(j+1,dm[1])  ]-v[i+dm[0]*WRAP(j-1,dm[1])  ])*td/2;
                j11     = (v[i+dm[0]*WRAP(j+1,dm[1])+m]-v[i+dm[0]*WRAP(j-1,dm[1])+m])*td/2 + 1.0;
                J0[o]   = j00*j11 - j10*j01;
            }
        }
        for(i=0; i<k; i++)
        {
            double *tmpp;
            composition_detjac(dm, t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    }
    else
    {
        for(j=0; j<dm[1]; j++)
        {
            for(i=0; i<dm[0]; i++)
            {
                int o   = i+dm[0]* j;
                t0[o  ] = (i+1) + v[o  ]*td;
                t0[o+m] = (j+1) + v[o+m]*td;
            }
        }
        for(i=0; i<k; i++)
        {
            double *tmpp;
            composition(dm, t0, t0, t1);
            tmpp = t0; t0   = t1; t1   = tmpp;
        }
    }

    if (optr != t0)
    {
        for(i=0; i<2*m; i++)
            t1[i] = t0[i];

        if (J0!=(double *)0)
            for(i=0; i<4*m; i++)
                J1[i] = J0[i];
    }
}

int pow2(int k)
{
    int j0, td = 1;
    for(j0=0; j0<k; j0++)
        td = td*2;
    return(td);
}

/*
 * J0 := J0*inv(I+diag(v0)*sc)
 */
void jac_div_smalldef(int dm[], double sc, double v0[], double J0[])
{
    int j0, j1;
    int m = dm[0]*dm[1];
    double sc2 = sc/2.0;
    double *v1 = v0+m;

    for(j1=0; j1<dm[1]; j1++)
    {
        int j1m1, j1p1;
        j1m1 = WRAP(j1-1,dm[1]);
        j1p1 = WRAP(j1+1,dm[1]);

        for(j0=0; j0<dm[0]; j0++)
        {
            int o, om1, op1;
            double j00,j01, j10,j11;
            double t00,t01, t10,t11;
            double idt;

            om1 = WRAP(j0-1,dm[0])+dm[0]*j1;
            op1 = WRAP(j0+1,dm[0])+dm[0]*j1;
            j00 = (v0[op1]-v0[om1])*sc2 + 1.0;
            j01 = (v1[op1]-v1[om1])*sc2;

            om1 = j0+dm[0]*j1m1;
            op1 = j0+dm[0]*j1p1;
            j10 = (v0[op1]-v0[om1])*sc2;
            j11 = (v1[op1]-v1[om1])*sc2 + 1.0;

            /*
            syms j00 j01 j10 j11
            syms t00 t01 t10 t11
            J1 = [j00 j01; j10 j11];
            J0 = [t00 t01; t10 t11];
            inv(J1)
            J1*J0
            */
            idt = 1.0/(j00*j11-j10*j01);
            t00 =  idt*j11;
            t01 = -idt*j10;
            t10 = -idt*j01;
            t11 =  idt*j00;


            o   = j0+dm[0]*j1;
            j00 = J0[o  ]; j01 = J0[o+m*2];
            j10 = J0[o+m]; j11 = J0[o+m*3];
            J0[o    ] = j00*t00+j01*t10;
            J0[o+m  ] = j10*t00+j11*t10;
            J0[o+m*2] = j00*t01+j01*t11;
            J0[o+m*3] = j10*t01+j11*t11;
        }
    }
}

/* Similar purpose to objfun, but uses a logistic regression model */
double initialise_objfun2(int dm[], double f[], double g[], double t0[], double J0[], double dj[], double b[], double A[])
{
    int j, m = dm[0]*dm[1];
    double ssl = 0.0, dt = 1.0;

    for(j=0; j<m; j++)
    {
        double x, y;
        int    ix, iy, ix1, iy1, k;
        double k11, k12, k21, k22;
        double dx0, dx1, dx2, dy0, dy1, dy2;
        double Y[128], dx[128], dy[128], sY;
        double ta11, ta12, ta22, tb1, tb2, tss;

        x    = t0[j  ]-1.0;
        y    = t0[j+m]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        ix   = WRAP(ix,dm[0]);
        iy   = WRAP(iy,dm[1]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);
        sY   = 0.0;
        for(k=0; k<dm[2]; k++)
        {
            k22   = f[ix  + dm[0]*iy +m*k];
            k12   = f[ix1 + dm[0]*iy +m*k];
            k21   = f[ix  + dm[0]*iy1+m*k];
            k11   = f[ix1 + dm[0]*iy1+m*k];

            Y[k]  = exp((k11*dx1 + k21*dx2)*dy1 + (k12*dx1 + k22*dx2)*dy2);
            dx0   =    ((k11     - k21    )*dy1 + (k12     - k22    )*dy2);
            dy0   =    ((k11*dx1 + k21*dx2)     - (k12*dx1 + k22*dx2)    );
            sY   += Y[k];
            dx[k] = J0[j    ]*dx0 + J0[j+  m]*dy0;
            dy[k] = J0[j+2*m]*dx0 + J0[j+3*m]*dy0;
        }
        ta11 = ta22 = ta12 = 0.0;
        tb1  = tb2  = 0.0;
        tss  = 0.0;
        for(k=0; k<dm[2]; k++)
        {
            double T = g[j + m*k], wt;
            int k1;
            Y[k] /= sY;
            tss  += log(Y[k])*T;
            tb1  += (Y[k]-T)*dx[k];
            tb2  += (Y[k]-T)*dy[k];
            for(k1=0; k1<k; k1++)
            {
                wt    = -Y[k]*Y[k1];
                ta11 += wt* dx[k]*dx[k1]*2.0;
                ta22 += wt* dy[k]*dy[k1]*2.0;
                ta12 += wt*(dx[k]*dy[k1]+dx[k1]*dy[k]);
            }
            wt    = Y[k]*(1.0-Y[k]);
            ta11 += wt*dx[k]*dx[k];
            ta22 += wt*dy[k]*dy[k];
            ta12 += wt*dx[k]*dy[k];
        }
        if (dj != (double *)0)
            dt = dj[j];

        A[j    ] = ta11*dt;
        A[j+  m] = ta22*dt;
        A[j+2*m] = ta12*dt;
        b[j  ]   = tb1*dt;
        b[j+m]   = tb2*dt;
        ssl     -= tss*dt;
    }
    return(ssl);
}

double initialise_objfun(int dm[], double f[], double g[], double t0[], double J0[], double dj[], double b[], double A[])
{
    int j, m = dm[0]*dm[1];
    double ssl = 0.0, dt = 1.0;

    if (dm[2]>1)
    {
        ssl = initialise_objfun2(dm, f, g, t0, J0, dj, b, A);
        return(ssl);
    }
    for(j=0; j<m; j++)
    {
        double x, y;
        int    ix, iy, ix1, iy1;
        double k11, k12, k21, k22;
        double dx0, dx1, dx2, dy0, dy1, dy2;
        double d, dx, dy;

        x    = t0[j  ]-1.0;
        y    = t0[j+m]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        ix   = WRAP(ix,dm[0]);
        iy   = WRAP(iy,dm[1]);
        ix1  = WRAP(ix+1,dm[0]);
        iy1  = WRAP(iy+1,dm[1]);

        k22  = f[ix  + dm[0]*iy ];
        k12  = f[ix1 + dm[0]*iy ];
        k21  = f[ix  + dm[0]*iy1];
        k11  = f[ix1 + dm[0]*iy1];

        d    = ((k11*dx1 + k21*dx2)*dy1 + (k12*dx1 + k22*dx2)*dy2) - g[j];
        dx0  = ((k11     - k21    )*dy1 + (k12     - k22    )*dy2);
        dy0  = ((k11*dx1 + k21*dx2)     - (k12*dx1 + k22*dx2)    );

        dx   = J0[j    ]*dx0 + J0[j+  m]*dy0;
        dy   = J0[j+2*m]*dx0 + J0[j+3*m]*dy0;

        if (dj != (double *)0)
            dt = dj[j];

        A[j    ] = dx*dx*dt;
        A[j+  m] = dy*dy*dt;
        A[j+2*m] = dx*dy*dt;

        b[j  ]   = dx*d*dt;
        b[j+m]   = dy*d*dt;

        ssl += d*d*dt;
    }
    return(0.5*ssl);
}


/*
 * t0 = Id + v0*sc
 * J0 = Id + I+diag(v0)*sc
 */
void smalldef_jac(int dm[], double sc, double v0[], double t0[], double J0[])
{
    int i, j;
    int m = dm[0]*dm[1];
    double sc2 = sc/2.0;
    double *v1 = v0+m;

    for(j=0; j<dm[1]; j++)
    {
        for(i=0; i<dm[0]; i++)
        {
            int o   = i+dm[0]*j;
            int om1, op1;
            t0[o  ] = (i+1) + v0[o]*sc;
            t0[o+m] = (j+1) + v1[o]*sc;

            om1 = WRAP(i-1,dm[0])+dm[0]*j;
            op1 = WRAP(i+1,dm[0])+dm[0]*j;
            J0[o    ] = (v0[op1]-v0[om1])*sc2 + 1.0;
            J0[o+  m] = (v1[op1]-v1[om1])*sc2;

            om1 = i+dm[0]*WRAP(j-1,dm[1]);
            op1 = i+dm[0]*WRAP(j+1,dm[1]);
            J0[o+2*m] = (v0[op1]-v0[om1])*sc2;
            J0[o+3*m] = (v1[op1]-v1[om1])*sc2 + 1.0;
        }
    }
}


void squaring(int dm[], int k, int save_transf, double b[], double A[], double t0[], double t1[], double J0[], double J1[])
{
    int i, j, m = dm[0]*dm[1];
    double *ptr = t0;

    for(i=0; i<k; i++)
    {
        double *buf1, *buf2;
        buf1 = t1; /* Re-use some memory */
        buf2 = J1;

        for(j=0; j<m; j++)
        {
            double tmp00, tmp01, tmp11, tmp1, tmp2, tmp3, tmp4;
            double x, y;
            double j11, j21, j12, j22, dt;

            x   = t0[j  ]-1.0;
            y   = t0[j+m]-1.0;
            j11 = J0[j  ]; j12 = J0[j+2*m];
            j21 = J0[j+m]; j22 = J0[j+3*m];
            dt  = j11*j22 - j12*j21;

         /* if (dt < 1e-9) dt = 1e-9;
            if (dt > 1e9 ) dt = 1e9; */

            tmp1      = samp(dm,b  ,x,y);
            tmp2      = samp(dm,b+m,x,y);

            buf1[j  ] = dt*(tmp1*j11+tmp2*j21);
            buf1[j+m] = dt*(tmp1*j12+tmp2*j22);

            tmp00 = samp(dm,A    ,x,y);
            tmp11 = samp(dm,A+  m,x,y);
            tmp01 = samp(dm,A+2*m,x,y);

            tmp1  = tmp00*j11+tmp01*j21;
            tmp2  = tmp01*j11+tmp11*j21;
            tmp3  = tmp00*j12+tmp01*j22;
            tmp4  = tmp01*j12+tmp11*j22;

         /* if (dt < 1e-9) dt = 1e-9; */
         /* if (dt > 1e9 ) dt = 1e9; */

            buf2[j    ] = dt*(tmp1*j11+tmp2*j21);
            buf2[j+  m] = dt*(tmp3*j12+tmp4*j22);
            buf2[j+2*m] = dt*(tmp1*j12+tmp2*j22);
        }
        for(j=0; j<2*m; j++) b[j] += buf1[j];
        for(j=0; j<3*m; j++) A[j] += buf2[j];

        if (save_transf || (i<k-1))
        {
            double *tmpp;
            composition_jacobian(dm, t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    }
    if (save_transf && ptr!=t0)
    {
        for(j=0; j<m*2; j++) t1[j] = t0[j];
        for(j=0; j<m*4; j++) J1[j] = J0[j];
    }
}

void unwrap(int dm[], double f[])
{
    int i,j;

    f[0] = f[0]-floor((f[0]-0)/dm[0]+0.5)*dm[0];
    for(j=0; j<dm[1]; j++)
    {
        double *pt = &f[j*dm[0]];
        if (j==0)
        {
            pt[0] = pt[0]-floor(pt[0]/dm[0]+0.5)*dm[0];
        }
        else
        {
            pt[0] = pt[0]-floor((pt[0]-pt[-dm[0]])/dm[0]+0.5)*dm[0];
        }
        for(i=1; i<dm[0]; i++)
        {
            pt[i] = pt[i]-floor((pt[i]-pt[i-1])/dm[0]+0.5)*dm[0];
        }
    }
    for(i=0; i<dm[0]; i++)
    {
        double *pt = &f[i+dm[0]*dm[1]];
        if (i==0)
        {
            pt[0] = pt[0]-floor(pt[0]/dm[1]+0.5)*dm[1];
        }
        else
        {
            pt[0] = pt[0]-floor((pt[0]-pt[-1])/dm[1]+0.5)*dm[1];
        }
        for(j=dm[0]; j<dm[1]*dm[0]; j+=dm[0])
        {
            pt[j] = pt[j]-floor((pt[j]-pt[j-dm[0]])/dm[1]+0.5)*dm[1];
        }
    }
}

int dartel_scratchsize(int dm[], int issym)
{
    int m1, m2;
    int m = dm[0]*dm[1];
    m1 = 15*m;
    if (issym) m1 += 5*m;
    m2 = 5*m+fmg2_scratchsize(dm);
    if (m1>m2)
        return(m1);
    else
        return(m2);
}

void dartel(int dm[], int k, double v[], double g[], double f[], double dj[], int rtype, double param[], double lmreg, int cycles, int nits, int issym, double ov[], double ll[], double *buf)
{
    double *sbuf;
    double *b, *A, *b1, *A1;
    double *t0, *t1, *J0, *J1;
    double sc;
    double ssl, ssp;
    double normb;
    int j, m = dm[0]*dm[1];

    /*
        Allocate memory.
          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
        [ A  A  A  t  t  J  J  J  J  t  t  J  J  J  J] for computing derivatives
        [ A  A  A s1 s2 s3 s4 s5 s6 s7 s8] for CGS solver
    */
    b    = ov;
    A    = buf;
    t0   = buf +  3*m;
    J0   = buf +  5*m;
    t1   = buf +  9*m;
    J1   = buf + 11*m;
    A1   = buf + 15*m;
    b1   = buf + 18*m;
    sbuf = buf +  3*m;

    sc = 1.0/pow2(k);

    expdef(dm, k, v, t0, t1, J0, J1);
    jac_div_smalldef(dm, sc, v, J0);
    ssl = initialise_objfun(dm, f, g, t0, J0, dj, b, A);
    smalldef_jac(dm, -sc, v, t0, J0);
    squaring(dm, k, issym, b, A, t0, t1, J0, J1);

    if (issym)
    {
        jac_div_smalldef(dm, -sc, v, J0);
        ssl += initialise_objfun(dm, g, f, t0, J0, (double *)0, b1, A1);
        smalldef_jac(dm, sc, v, t0, J0);
        squaring(dm, k, 0, b1, A1, t0, t1, J0, J1);
        for(j=0; j<m*2; j++) b[j] -= b1[j];
        for(j=0; j<m*3; j++) A[j] += A1[j];
    }

    if (rtype==0)
        LtLf_le(dm, v, param, t1);
    else if (rtype==1)
        LtLf_me(dm, v, param, t1);
    else
        LtLf_be(dm, v, param, t1);

    ssp = 0.0;
    for(j=0; j<2*m; j++)
    {
        b[j] = b[j]*sc + t1[j];
        ssp += t1[j]*v[j];
    }
    normb = norm(2*m,b);

    for(j=0; j<3*m; j++) A[j] *= sc;
    for(j=0; j<2*m; j++) A[j] += lmreg;

    /* Solve equations for Levenberg-Marquardt update:
     * v = v - inv(H + L'*L + R)*(d + L'*L*v)
     *     v: velocity or flow field
     *     H: matrix of second derivatives
     *     L: regularisation (L'*L is the inverse of the prior covariance)
     *     R: Levenberg-Marquardt regularisation
     *     d: vector of first derivatives
     */
 /* cgs2(dm, A, b, rtype, param, 1e-8, 4000, sbuf, sbuf+2*m, sbuf+4*m, sbuf+6*m); */
    fmg2(dm, A, b, rtype, param, cycles, nits, sbuf, sbuf+2*m);

    for(j=0; j<2*m; j++) ov[j] = v[j] - sbuf[j];
    ll[0] = ssl;
    ll[1] = ssp*0.5;
    ll[2] = normb;
}

