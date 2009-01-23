/* $Id: optimN.c 2644 2009-01-23 13:01:50Z john $ */
/* (c) John Ashburner (2007) */

#include<mex.h>
#include<math.h>
extern double log(double x);
#define MAXD3 128

#include "optimN.h"

#ifdef NEUMANN
    /* Neumann boundary condition */
    static int neumann(int i,  int m)
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

static void choldc(int n, double a[], double p[])
{
    int i, j, k;
    double sm, sm0;

    sm0  = 1e-16;
    for(i=0; i<n; i++) sm0 = sm0 + a[i*n+i];
    sm0 *= 1e-4;
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

static void Atimesp1(int dm[], float A[], float p[], float Ap[])
{
    int i,j, m = dm[0]*dm[1]*dm[2];
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
        int k, o;
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

static void get_a(int dm3, int i, float *pa[],  double a[])
{
    int m, n;
    int o = dm3;
    for(m=0; m<dm3; m++)
    {
        a[m+dm3*m] = pa[m][i];
        for(n=m+1; n<dm3; n++,o++)
            a[n+dm3*m] = a[m+dm3*n] = pa[o][i];
    }
}

static double sumsq_me(int dm[], float a[], float b[], double s[], double scal[], float u[])
{
    double w000, w001, w010, w100;
    double ss = 0.0;
    int i, j, k, m;

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
            float *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
            double a1[MAXD3*MAXD3];
            int jm1,jp1,im1,ip1;

            for(m=0; m<dm[3]; m++)
            {
                pu[m]  = u+dm[0]*(j+dm[1]*(k+dm[2]*m));
                pb[m]  = b+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }
            for(m=0; m<(dm[3]*(dm[3]+1))/2; m++)
                pa[m]  = a+dm[0]*(j+dm[1]*(k+dm[2]*m));

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                float *p[MAXD3];
                double tmp;

                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;

                for(m=0; m<dm[3]; m++) p[m] = &(pu[m][i]);
                get_a(dm[3], i, pa, a1);
                for(m=0; m<dm[3]; m++)
                {
                    int n;
                    float *pm = &pu[m][i];
                    double *a11 = a1 + dm[3]*m;
                    tmp = (w000* pm[0] + 
                         + w001*(pm[km1] + pm[kp1])
                         + w010*(pm[jm1] + pm[jp1])
                         + w100*(pm[im1] + pm[ip1]))*scal[m]
                         - pb[m][i];
                    for(n=0; n<dm[3]; n++)
                        tmp += a11[n]*pu[n][i];
                    ss += tmp*tmp;
                }
            }
        }
    }
    return(ss);
}

void LtLf_me(int dm[], float f[], double s[], double scal[], float g[])
{
    double w000,w001,w010,w100;
    int k;

    w000 = s[3]*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + s[5];
    w001 = s[3]*(-s[2]*s[2]);
    w010 = s[3]*(-s[1]*s[1]);
    w100 = s[3]*(-s[0]*s[0]);

    for(k=0; k<dm[2]; k++)
    {
        int km1, kp1, j;
        km1 = (BOUND(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (BOUND(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            float *pf[MAXD3], *pg[MAXD3];
            int jm1,jp1,im1,ip1, m;

            for(m=0; m<dm[3]; m++)
            {
                pf[m]  = f+dm[0]*(j+dm[1]*(k+dm[2]*m));
                pg[m]  = g+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }

            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

            for(m=0; m<dm[3]; m++)
            {
                int i;
                float *pf1 = pf[m], *pg1 = pg[m];
                for(i=0; i<dm[0]; i++)
                {
                    float *p = pf1 + i;
                    im1    = BOUND(i-1,dm[0])-i;
                    ip1    = BOUND(i+1,dm[0])-i;
                    pg1[i] = (w000*p[0] + w001*(p[km1] + p[kp1]) + w010*(p[jm1] + p[jp1]) + w100*(p[im1] + p[ip1]))*scal[m];
                }
            }
        }
    }
}

static void relax_me(int dm[], float a[], float b[], double s[], double scal[], int nit, float u[])
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
        int m;

#       ifdef VERBOSE
            printf(" %g", sumsq_me(dm, a, b, s, scal, u));
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
                float *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
                double a1[MAXD3*MAXD3], cp[MAXD3], su[MAXD3];
                int jm1,jp1,im1,ip1;

                for(m=0; m<dm[3]; m++)
                {
                    pu[m]  = u+dm[0]*(j+dm[1]*(k+dm[2]*m));
                    pb[m]  = b+dm[0]*(j+dm[1]*(k+dm[2]*m));
                }
                for(m=0; m<(dm[3]*(dm[3]+1))/2; m++)
                    pa[m]  = a+dm[0]*(j+dm[1]*(k+dm[2]*m));

                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];

                istart = (jstart == (j%2));

                for(i=istart; i<dm[0]; i+=2)
                {
                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;

                    get_a(dm[3], i, pa, a1);
                    for(m=0; m<dm[3]; m++)
                    {
                        float *pm = &pu[m][i];
                        su[m]     = pb[m][i]-(w001*(pm[km1] + pm[kp1]) + w010*(pm[jm1] + pm[jp1]) + w100*(pm[im1] + pm[ip1]))*scal[m];
                        a1[m+dm[3]*m] += w000*scal[m];
                    }
                    choldc(dm[3],a1,cp);
                    cholls(dm[3],a1,cp,su,su);
                    for(m=0; m<dm[3]; m++) pu[m][i] = su[m];
                }
            }
        }
    }
#   ifdef VERBOSE
        printf(" %g\n", sumsq_me(dm, a, b, s, scal, u));
#   endif
}

static void Atimesp_me(int dm[], float A[], double param[], double scal[], float p[], float Ap[])
{
    LtLf_me(dm, p, param, scal, Ap);
    Atimesp1(dm, A, p, Ap);
}

static double sumsq_be(int dm[], float a[], float b[], double s[], double scal[], float u[])
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
            int i,m, jm2,jm1,jp1,jp2;
            float *p[MAXD3], *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
            double a1[MAXD3*MAXD3];

            for(m=0; m<dm[3]; m++)
            {
                pu[m]  = u+dm[0]*(j+dm[1]*(k+dm[2]*m));
                pb[m]  = b+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }
            for(m=0; m<(dm[3]*(dm[3]+1))/2; m++)
                pa[m]  = a+dm[0]*(j+dm[1]*(k+dm[2]*m));

            jm2 = (BOUND(j-2,dm[1])-j)*dm[0];
            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];
            jp2 = (BOUND(j+2,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                int m, im2,im1,ip1,ip2;
                double tmp;

                im2 = BOUND(i-2,dm[0])-i;
                im1 = BOUND(i-1,dm[0])-i;
                ip1 = BOUND(i+1,dm[0])-i;
                ip2 = BOUND(i+2,dm[0])-i;

                for(m=0; m<dm[3]; m++) p[m] = &(pu[m][i]);
                get_a(dm[3], i, pa, a1);
                for(m=0; m<dm[3]; m++)
                {
                    int n;
                    float *pm = p[m];
                    double *a11 = a1 + dm[3]*m;
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

                    for(n=0; n<dm[3]; n++) tmp += a11[n]*p[n][0];
                    ss += tmp*tmp;
                }
            }
        }
    }
    return(ss);
}

void LtLf_be(int dm[], float f[], double s[], double scal[], float g[])
{
    int k;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    /*
        syms s0 s1 s2 s3 s4
        zz = sym(zeros(3,3));
        K1 = cat(3,zz,[0 -s0*s0 0; 0 2*s0*s0 0; 0 -s0*s0 0],zz);
        K2 = cat(3,zz,[0 0 0; -s1*s1 2*s1*s1 -s1*s1; 0 0 0],zz);
        K3 = sym(zeros(3,3,3));
        K3(2,2,1) = -s2*s2;
        K3(2,2,2) = 2*s2*s2;
        K3(2,2,3) = -s2*s2;

        K  = K1+K2+K3;
        K(2,2,2) = K(2,2,2) + s4;
        % L  = convn(K,K)
        L  = sym(zeros(5,5,5));
        for i=1:3,
            for j=1:3,
                for k=1:3,
                    L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) = L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) + K(i,j,k)*K;
                end;
            end;
        end;
        disp(s3*L(3:end,3:end,3:end))
    */

    w000 = s[3]*(6*(s[0]*s[0]*s[0]*s[0]+s[1]*s[1]*s[1]*s[1]+s[2]*s[2]*s[2]*s[2])
                +8*(s[0]*s[0]*s[1]*s[1]+s[0]*s[0]*s[2]*s[2]+s[1]*s[1]*s[2]*s[2])
                +4*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])*s[4] + s[4]*s[4]) + s[5];
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
            int i,m, jm2,jm1,jp1,jp2;
            float *pf[MAXD3], *pg[MAXD3];

            for(m=0; m<dm[3]; m++)
            {
                pf[m]  = f+dm[0]*(j+dm[1]*(k+dm[2]*m));
                pg[m]  = g+dm[0]*(j+dm[1]*(k+dm[2]*m));
            }

            jm2 = (BOUND(j-2,dm[1])-j)*dm[0];
            jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
            jp1 = (BOUND(j+1,dm[1])-j)*dm[0];
            jp2 = (BOUND(j+2,dm[1])-j)*dm[0];

            for(m=0; m<dm[3]; m++)
            {
                int im2,im1,ip1,ip2;
                float *pf1 = pf[m], *pg1 = pg[m];
                for(i=0; i<dm[0]; i++)
                {
                    float *p = &pf1[i];
                    im2 = BOUND(i-2,dm[0])-i;
                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;
                    ip2 = BOUND(i+2,dm[0])-i;

                    pg1[i] =(w000* p[0]
                           + w010*(p[    jm1    ] + p[    jp1    ])
                           + w020*(p[    jm2    ] + p[    jp2    ])
                           + w100*(p[im1        ] + p[ip1        ])
                           + w110*(p[im1+jm1    ] + p[ip1+jm1    ] + p[im1+jp1    ] + p[ip1+jp1    ])
                           + w200*(p[im2        ] + p[ip2        ])
                           + w001*(p[        km1] + p[        kp1])
                           + w101*(p[im1    +km1] + p[ip1    +km1] + p[im1    +kp1] + p[ip1    +kp1])
                           + w011*(p[    jm1+km1] + p[    jp1+km1] + p[    jm1+kp1] + p[    jp1+kp1])
                           + w002*(p[        km2] + p[        kp2]))*scal[m];
                }
            }
        }
    }
}

static void relax_be(int dm[], float a[], float b[], double s[], double scal[], int nit, float u[])
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
       This is an attempt to stabilise the relaxation. */
    reg = 1.0001*(2.0*(w200+w020+w002)-2.0*(w100+w010+w001)+4.0*(w110+w011+w101)) - w000;
    if (reg<0.0) reg = w000*0.00001;

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("%dx%dx%d (%g): ", dm[0],dm[1],dm[2],reg);
        printf(" %g", sumsq_be(dm, a, b, s, scal, u));
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
                float *pu[MAXD3], *pb[MAXD3], *pa[(MAXD3*(MAXD3+1))/2];
                double a1[MAXD3*MAXD3], cp[MAXD3], su[MAXD3];
                int m, jm2,jm1,jp1,jp2;

                for(m=0; m<dm[3]; m++)
                {
                    pu[m]  = u+dm[0]*(j+dm[1]*(k+dm[2]*m));
                    pb[m]  = b+dm[0]*(j+dm[1]*(k+dm[2]*m));
                }
                for(m=0; m<(dm[3]*(dm[3]+1))/2; m++)
                    pa[m]  = a+dm[0]*(j+dm[1]*(k+dm[2]*m));

                jm2 = (BOUND(j-2,dm[1])-j)*dm[0];
                jm1 = (BOUND(j-1,dm[1])-j)*dm[0];
                jp1 = (BOUND(j+1,dm[1])-j)*dm[0];
                jp2 = (BOUND(j+2,dm[1])-j)*dm[0];

                for(i=istart; i!=iend; i+=iskip)
                {
                    int im2,im1,ip1,ip2;

                    im2 = BOUND(i-2,dm[0])-i;
                    im1 = BOUND(i-1,dm[0])-i;
                    ip1 = BOUND(i+1,dm[0])-i;
                    ip2 = BOUND(i+2,dm[0])-i;

                    get_a(dm[3], i, pa, a1);
                    for(m=0; m<dm[3]; m++)
                    {
                        int n;
                        float *pm = &pu[m][i];
                        su[m]     = pb[m][i]-
                                  ( w000* pm[0]
                                  + w010*(pm[    jm1    ] + pm[    jp1    ])
                                  + w020*(pm[    jm2    ] + pm[    jp2    ])
                                  + w100*(pm[im1        ] + pm[ip1        ])
                                  + w110*(pm[im1+jm1    ] + pm[ip1+jm1    ] + pm[im1+jp1    ] + pm[ip1+jp1    ])
                                  + w200*(pm[im2        ] + pm[ip2        ])
                                  + w001*(pm[        km1] + pm[        kp1])
                                  + w101*(pm[im1    +km1] + pm[ip1    +km1] + pm[im1    +kp1] + pm[ip1    +kp1])
                                  + w011*(pm[    jm1+km1] + pm[    jp1+km1] + pm[    jm1+kp1] + pm[    jp1+kp1])
                                  + w002*(pm[        km2] + pm[        kp2]))*scal[m];
                        for(n=0; n<dm[3]; n++) su[m] -= a1[m*dm[3]+n]*pu[n][i];
                        a1[m+dm[3]*m] += (w000+reg)*scal[m];
                    }
                    choldc(dm[3],a1,cp);
                    cholls(dm[3],a1,cp,su,su);
                    for(m=0; m<dm[3]; m++) pu[m][i] += su[m];

                    /* ss  += sux*sux + suy*suy + suz*suz; */
                }
            }
        }
#       ifdef VERBOSE
            printf(" %g", sumsq_be(dm, a, b, s, scal, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}


static void Atimesp_be(int dm[], float A[], double param[], double scal[], float p[], float Ap[])
{
    LtLf_be(dm, p, param, scal, Ap);
    Atimesp1(dm, A, p, Ap);
}

float norm(int m, float a[])
{
    int i;
    double dp = 0.0;
    for(i=0; i<m; i++)
        dp += a[i]*a[i];
    return(sqrt(dp));
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

int fmg_scratchsize(int n0[])
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
    return((n0[3]*n0[0]*n0[1]*n0[2] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + (n0[3]*3+(n0[3]*(n0[3]+1))/2)*bs));
}

/*
    Full Multigrid solver.  See Numerical Recipes (second edition) for more
    information
*/
void fmg(int n0[], float *a0, float *b0, int rtype, double param0[], double scal[], int c, int nit,
          float *u0, float *scratch)
{
    int i, j, ng, bs;
     int n[32][4], m[32];
    float *bo[32], *a[32], *b[32], *u[32], *res, *rbuf;
    double param[32][6];
    void (*relax)(), (*Atimesp)();
    double (*sumsq)();

    if (rtype == 1)
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
        restrict(n0[3],n[j-1],bo[j-1],n[j],bo[j],rbuf);
        restrict((n0[3]*(n0[3]+1))/2,n[j-1],a[j-1],n[j],a[j],rbuf);

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

                restrict(n0[3],n[jj],res,n[jj+1],b[jj+1],rbuf);
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

/*  printf("end=%g\n", sumsq(n0, a0, b0, param[0], u0)); */
}

