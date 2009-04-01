/* $Id: optim1.c 3030 2009-04-01 13:51:00Z guillaume $ */
/* (c) John Ashburner (2007) */

#include<math.h>
#include "mex.h"
#include "optim1.h"
#define WRAP(i,m) (((i)>=0) ? (i)%(m) : ((m)+(i)%(m))%m)
extern double sqrt(double x), log(double x);

static void Atimesp1(int dm[], double A[], double p[], double Ap[])
{
    int i, m = dm[0]*dm[1]*dm[2];
    for(i=0; i<m; i++)
        Ap[i] += A[i]*p[i];
}

double sumsq_me(int dm[], double a[], double b[], double s[], double u[])
{
    double w000, w001, w010, w100;
    double ss = 0.0;
    int i, j, k;

    w000 = s[3]*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + s[4];
    w001 = s[3]*(-s[2]*s[2]);
    w010 = s[3]*(-s[1]*s[1]);
    w100 = s[3]*(-s[0]*s[0]);

    for(k=0; k<dm[2]; k++)
    {
        int km1,kp1;
        km1 = (WRAP(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (WRAP(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            double *pux, *pbx, *paxx;
            int jm1,jp1,im1,ip1;

            pux  = u+dm[0]*(j+dm[1]*k);
            pbx  = b+dm[0]*(j+dm[1]*k);
            paxx = a+dm[0]*(j+dm[1]*k);

            jm1 = (WRAP(j-1,dm[1])-j)*dm[0];
            jp1 = (WRAP(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                double *px = &pux[i];
                double tmp;

                im1 = WRAP(i-1,dm[0])-i;
                ip1 = WRAP(i+1,dm[0])-i;

                tmp = (w000+paxx[i])*px[0]
                     + w001*(px[km1] + px[kp1])
                     + w010*(px[jm1] + px[jp1])
                     + w100*(px[im1] + px[ip1])
                     - pbx[i];
                ss += tmp*tmp;
            }
        }
    }
    return(ss);
}

void LtLf_me(int dm[], double f[], double s[], double g[])
{
    int i, j, k, km1,kp1, jm1,jp1, im1,ip1;
    double *pgx, *pfx;
    double w000,w001,w010,w100;

    w000 = s[3]*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + s[4];
    w001 = s[3]*(-s[2]*s[2]);
    w010 = s[3]*(-s[1]*s[1]);
    w100 = s[3]*(-s[0]*s[0]);

    for(k=0; k<dm[2]; k++)
    {
        km1 = (WRAP(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (WRAP(k+1,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            pgx = g+dm[0]*(j+dm[1]*k);
            pfx = f+dm[0]*(j+dm[1]*k);
            jm1 = (WRAP(j-1,dm[1])-j)*dm[0];
            jp1 = (WRAP(j+1,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                double *px = &pfx[i];
                im1 = WRAP(i-1,dm[0])-i;
                ip1 = WRAP(i+1,dm[0])-i;
                pgx[i] = w000*px[0] + w001*(px[km1] + px[kp1]) + w010*(px[jm1] + px[jp1]) + w100*(px[im1] + px[ip1]);
            }
        }
    }
}

static void relax_me(int dm[], double a[], double b[], double s[], int nit, double u[])
{
    int it;
    double w000,w001,w010,w100;

    w000 = s[3]*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + s[4];
    w001 = s[3]*(-s[2]*s[2]);
    w010 = s[3]*(-s[1]*s[1]);
    w100 = s[3]*(-s[0]*s[0]);

#ifdef VERBOSE
    for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
    printf("%dx%dx%d: ", dm[0],dm[1],dm[2]);
#endif

    for(it=0; it<2*nit; it++)
    {
        int k, kstart;
        int j, jstart;
        int i, istart;

#ifdef VERBOSE
        printf(" %g", sumsq_me(dm, a, b, s, u));
#endif

        kstart = it%2;
        for(k=0; k<dm[2]; k++)
        {
            int km1, kp1;
            km1 = (WRAP(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (WRAP(k+1,dm[2])-k)*dm[0]*dm[1];

            jstart = (kstart == (k%2));
            for(j=0; j<dm[1]; j++)
            {
                double *pux, *pbx, *paxx;
                int jm1,jp1, im1,ip1;

                pux  = u+dm[0]*(j+dm[1]*k);
                pbx  = b+dm[0]*(j+dm[1]*k);
                paxx = a+dm[0]*(j+dm[1]*k);

                jm1 = (WRAP(j-1,dm[1])-j)*dm[0];
                jp1 = (WRAP(j+1,dm[1])-j)*dm[0];

                istart = (jstart == (j%2));

                for(i=istart; i<dm[0]; i+=2)
                {
                    double sux;
                    double *px = pux+i;
                    im1 = WRAP(i-1,dm[0])-i;
                    ip1 = WRAP(i+1,dm[0])-i;
                    sux = pbx[i]-(w001*(px[km1] + px[kp1]) + w010*(px[jm1] + px[jp1]) + w100*(px[im1] + px[ip1]));
                    *px = sux/(paxx[i]+w000);
                }
            }
        }
    }
#ifdef VERBOSE
    printf(" %g\n", sumsq_me(dm, a, b, s, u));
#endif
}

static void Atimesp_me(int dm[], double A[], double param[], double p[], double Ap[])
{
    LtLf_me(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

double sumsq_be(int dm[], double a[], double b[], double s[], double u[])
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
                +8*(s[0]*s[0]*s[1]*s[1]+s[0]*s[0]*s[2]*s[2]+s[1]*s[1]*s[2]*s[2])) + s[4];
    w100 = s[3]*(-4*s[0]*s[0]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w200 = s[3]*s[0]*s[0]*s[0]*s[0];
    w010 = s[3]*(-4*s[1]*s[1]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w110 = s[3]*2*s[0]*s[0]*s[1]*s[1];
    w020 = s[3]*s[1]*s[1]*s[1]*s[1];

    w001 = s[3]*(-4*s[2]*s[2]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w101 = s[3]*2*s[0]*s[0]*s[2]*s[2];
    w011 = s[3]*2*s[1]*s[1]*s[2]*s[2];

    w002 = s[3]*s[2]*s[2]*s[2]*s[2];

    for(k=0; k<dm[2]; k++)
    {
        int j, km2,km1,kp1,kp2;
        km2 = (WRAP(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (WRAP(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (WRAP(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (WRAP(k+2,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            double *pux, *pbx, *paxx;
            int i, jm2,jm1,jp1,jp2;

            pux  = u+dm[0]*(j+dm[1]*k);
            pbx  = b+dm[0]*(j+dm[1]*k);
            paxx = a+dm[0]*(j+dm[1]*k);
            jm2  = (WRAP(j-2,dm[1])-j)*dm[0];
            jm1  = (WRAP(j-1,dm[1])-j)*dm[0];
            jp1  = (WRAP(j+1,dm[1])-j)*dm[0];
            jp2  = (WRAP(j+2,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                int im2,im1,ip1,ip2;
                double *px = pux+i;
                double tmp;

                im2 = WRAP(i-2,dm[0])-i;
                im1 = WRAP(i-1,dm[0])-i;
                ip1 = WRAP(i+1,dm[0])-i;
                ip2 = WRAP(i+2,dm[0])-i;

                tmp = (w000+paxx[i])*px[0]
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
            }
        }
    }
    return(ss);
}

void LtLf_be(int dm[], double f[], double s[], double g[])
{
    int k;
    double w000,w100,w200,
          w010,w110,
          w020,
          w001,w101,
          w011,
          w002;

    w000 = s[3]*(6*(s[0]*s[0]*s[0]*s[0]+s[1]*s[1]*s[1]*s[1]+s[2]*s[2]*s[2]*s[2])
                +8*(s[0]*s[0]*s[1]*s[1]+s[0]*s[0]*s[2]*s[2]+s[1]*s[1]*s[2]*s[2])) + s[4];
    w100 = s[3]*(-4*s[0]*s[0]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w200 = s[3]*s[0]*s[0]*s[0]*s[0];
    w010 = s[3]*(-4*s[1]*s[1]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w110 = s[3]*2*s[0]*s[0]*s[1]*s[1];
    w020 = s[3]*s[1]*s[1]*s[1]*s[1];

    w001 = s[3]*(-4*s[2]*s[2]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w101 = s[3]*2*s[0]*s[0]*s[2]*s[2];
    w011 = s[3]*2*s[1]*s[1]*s[2]*s[2];

    w002 = s[3]*s[2]*s[2]*s[2]*s[2];

    for(k=0; k<dm[2]; k++)
    {
        int j, km2,km1,kp1,kp2;
        km2 = (WRAP(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (WRAP(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (WRAP(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (WRAP(k+2,dm[2])-k)*dm[0]*dm[1];

        for(j=0; j<dm[1]; j++)
        {
            int i, jm2,jm1,jp1,jp2;
            double *pgx, *pfx;

            pgx = g+dm[0]*(j+dm[1]*k);
            pfx = f+dm[0]*(j+dm[1]*k);
            jm2 = (WRAP(j-2,dm[1])-j)*dm[0];
            jm1 = (WRAP(j-1,dm[1])-j)*dm[0];
            jp1 = (WRAP(j+1,dm[1])-j)*dm[0];
            jp2 = (WRAP(j+2,dm[1])-j)*dm[0];

            for(i=0; i<dm[0]; i++)
            {
                int im2,im1,ip1,ip2;
                double *px = &pfx[i];

                im2 = WRAP(i-2,dm[0])-i;
                im1 = WRAP(i-1,dm[0])-i;
                ip1 = WRAP(i+1,dm[0])-i;
                ip2 = WRAP(i+2,dm[0])-i;

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
            }
        }
    }
}

static void relax_be(int dm[], double a[], double b[], double s[], int nit, double u[])
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
                +8*(s[0]*s[0]*s[1]*s[1]+s[0]*s[0]*s[2]*s[2]+s[1]*s[1]*s[2]*s[2])) + s[4];
    w100 = s[3]*(-4*s[0]*s[0]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w200 = s[3]*s[0]*s[0]*s[0]*s[0];
    w010 = s[3]*(-4*s[1]*s[1]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w110 = s[3]*2*s[0]*s[0]*s[1]*s[1];
    w020 = s[3]*s[1]*s[1]*s[1]*s[1];

    w001 = s[3]*(-4*s[2]*s[2]*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
    w101 = s[3]*2*s[0]*s[0]*s[2]*s[2];
    w011 = s[3]*2*s[1]*s[1]*s[2]*s[2];

    w002 = s[3]*s[2]*s[2]*s[2]*s[2];

    /* For stability in Gauss-Seidel relaxation, the magnitude of the diagonal element must
       exceed the sum of the magnitudes of the off diagonal elements of each column or row
       (see e.g. http://www.mathpages.com/home/kmath175/kmath175.htm). */
    reg = (2.0*(w200+w020+w002)-2.0*(w100+w010+w001)+4.0*(w110+w011+w101)) - w000;
    if (reg<0.0) reg = 0.0;
#ifdef VERBOSE
    for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
    printf("%dx%dx%d (%g): ", dm[0],dm[1],dm[2],reg);
#endif

    for(it=0; it<8*nit; it++)
    {
        int k, kstart,kend,kskip;
        int j, jstart,jend,jskip;
        int i, istart,iend,iskip;
        double ss = 0.0;

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
            km2 = (WRAP(k-2,dm[2])-k)*dm[0]*dm[1];
            km1 = (WRAP(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (WRAP(k+1,dm[2])-k)*dm[0]*dm[1];
            kp2 = (WRAP(k+2,dm[2])-k)*dm[0]*dm[1];

            for(j=jstart; j!=jend; j+=jskip)
            {
                double *pux, *pbx, *paxx;
                int jm2,jm1,jp1,jp2;

                pux  = u+dm[0]*(j+dm[1]* k);
                pbx  = b+dm[0]*(j+dm[1]* k);
                paxx = a+dm[0]*(j+dm[1]* k);
                jm2  = (WRAP(j-2,dm[1])-j)*dm[0];
                jm1  = (WRAP(j-1,dm[1])-j)*dm[0];
                jp1  = (WRAP(j+1,dm[1])-j)*dm[0];
                jp2  = (WRAP(j+2,dm[1])-j)*dm[0];

                for(i=istart; i!=iend; i+=iskip)
                {
                    int im2,im1,ip1,ip2;
                    double sux;
                    double *px = pux+i;

                    im2 = WRAP(i-2,dm[0])-i;
                    im1 = WRAP(i-1,dm[0])-i;
                    ip1 = WRAP(i+1,dm[0])-i;
                    ip2 = WRAP(i+2,dm[0])-i;

                    sux = pbx[i] -((w000+paxx[i])*px[0]
                                  + w010*(px[    jm1    ] + px[    jp1    ])
                                  + w020*(px[    jm2    ] + px[    jp2    ])
                                  + w100*(px[im1        ] + px[ip1        ])
                                  + w110*(px[im1+jm1    ] + px[ip1+jm1    ] + px[im1+jp1    ] + px[ip1+jp1    ])
                                  + w200*(px[im2        ] + px[ip2        ])
                                  + w001*(px[        km1] + px[        kp1])
                                  + w101*(px[im1    +km1] + px[ip1    +km1] + px[im1    +kp1] + px[ip1    +kp1])
                                  + w011*(px[    jm1+km1] + px[    jp1+km1] + px[    jm1+kp1] + px[    jp1+kp1])
                                  + w002*(px[        km2] + px[        kp2]));
                    ss  += sux*sux;
                    *px += sux/(paxx[i] + w000 + reg);
                }
            }
        }
#ifdef VERBOSE
        printf(" %g", sumsq_be(dm, a, b, s, u));
#endif
    }
#ifdef VERBOSE
    printf("\n"); 
#endif
}


static void Atimesp_be(int dm[], double A[], double param[], double p[], double Ap[])
{
    LtLf_be(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

static double dotprod(int m, double a[], double b[])
{
    int i;
    double dp = 0.0;
    for(i=0; i<m; i++)
        dp += a[i]*b[i];
    return(dp);
}

static void addscaled(int m, double a[], double b[], double s)
{
    int i;
    for(i=0; i<m; i++)
        a[i] += s*b[i];
}

double norm(int m, double a[])
{
    return(sqrt(dotprod(m, a, a)));
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

void cgs(int dm[], double A[], double b[], int rtype, double param[], double tol, int nit,
             double x[], double r[], double p[], double Ap[])
{
    int i, m = dm[0]*dm[1]*dm[2], it;
    double rtr, nb, rtrold, alpha, beta;
    void (*Atimesp)();

    /* printf("\n **** %dx%d ****\n",dm[0],dm[1]); */
    if (rtype == 1)
        Atimesp = Atimesp_me;
    else
        Atimesp = Atimesp_be;

    nb      = tol*norm(m,b);

    if (0)
    {
        /* Assuming starting estimates of zeros */
        /* x    = zeros(size(b)); */
        for(i=0; i<m;i++)
            x[i] = 0.0;

        /* r    = b; */
        for(i=0; i<m;i++)
            r[i] = b[i];
    }
    else
    {
        /* Assume starting estimates are passed as arguments */
        /* r    = b-A*x; */
        Atimesp(dm, A, param, x, Ap);
        for(i=0; i<m;i++)
            r[i] = b[i]-Ap[i];
    }

    /* rtr  = r'*r; */
    rtr     = dotprod(m, r, r);

    /* p    = zeros(size(b)); */
    for(i=0; i<m;i++) p[i] = 0.0;

    /* beta = 0; */
    beta    = 0.0;

    /* for it=1:nit, */
    for(it=0; it<nit; it++)
    {
        /* if norm(r) < tol*norm(b), break; end; */
        if (norm(m,r) < nb) break;

        /* p      = r + beta*p; */
        for(i=0; i<m; i++) p[i]  = r[i] + beta*p[i];

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

static void resized_plane(int na[], double *a, int nc[], double *c, double *b)
{
    int i, j, o,om,op;
    double loc, s, w, wm, wp;
    double *ap, *bp, *cp;
    /* a - na[0]*na[1]
     * c - nc[0]*nc[1]
     * b - na[0]*nc[1]
     */

    s = (double)na[1]/(double)nc[1];
    for(j=0; j<nc[1]; j++)
    {
        loc = (j+0.5)*s-0.5;
        o   = floor(loc+0.5);
        om  = WRAP(o-1,na[1])*na[0];
        op  = WRAP(o+1,na[1])*na[0];
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
        om  = WRAP(o-1,na[0]);
        op  = WRAP(o+1,na[0]);
        w   = wt2( o   -loc);
        wp  = wt2((o+1)-loc);
        wm  = wt2((o-1)-loc);
        for(bp=b, cp=c+i, ap=bp+na[0]*nc[1]; bp<ap; bp+=na[0], cp+=nc[0])
            *cp = wm*bp[om]+w*bp[o]+wp*bp[op];
    }
}

void resize(int na[], double *a, int nc[], double *c, double *b)
{
    int j, k, o=-999999,om,op, m, oo;
    double loc, s, w, wm, wp;
    double *bp, *cp, *pl[3];
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
        om     = WRAP(o-1,na[2]);
        op     = WRAP(o+1,na[2]);

        if (o==oo)
        {   /* do nothing */
        }
        else if (o==oo+1)
        {   /* Shift by 1 */
            double *tp;
            tp    = pl[0];
            pl[0] = pl[1];
            pl[1] = pl[2];
            pl[2] = tp;
            resized_plane(na, a+na[0]*na[1]*op,nc,pl[2],bp);
        }
        else if (o==oo+2)
        {   /* Shift by 2 */
            double *tp;
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

static void zeros(int n, double *a)
{
    int i;
    for(i=0; i<n; i++)
        a[i] = 0.0;
}

static void copy(int n, double *a, double *b)
{
    int i;
    for(i=0; i<n; i++)
        b[i] = a[i];
}

static void addto(int n, double *a, double *b)
{
    int i;
    for(i=0; i<n; i++)
        a[i] += b[i];
}

int fmg_scratchsize(int n0[])
{
    int    n[32][3], m[32], bs=0, j;
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
    return((n0[0]*n0[1]*n0[2] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 4*bs));
}

/*
    Full Multigrid solver.  See Numerical Recipes (second edition) for more information
*/
void fmg(int n0[], double *a0, double *b0, int rtype, double param0[], int c, int nit,
          double *u0, double *scratch)
{
    int i, j, ng, bs;
    int    n[32][3], m[32];
    double *bo[32], *a[32], *b[32], *u[32], *res, *rbuf, param[32][5];
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

#ifdef VERBOSE
    printf("start=%g\n", sumsq(n0, a0, b0, param0, u0));
#endif

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
    rbuf   = scratch + m[0];
    bo[1]  = scratch + m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 0*bs;
    b[1]   = scratch + m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 1*bs;
    u[1]   = scratch + m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 2*bs;
    a[1]   = scratch + m[0] + n[0][0]*n[1][1]+3*n[0][0]*n[0][1] + 3*bs;

    for(j=2; j<ng; j++)
    {
        bo[j] = bo[j-1]+m[j-1];
        b[j]  =  b[j-1]+m[j-1];
        u[j]  =  u[j-1]+m[j-1];
        a[j]  =  a[j-1]+m[j-1];
    }

    for(j=1; j<ng; j++)
    {
        resize(n[j-1],bo[j-1],n[j],bo[j],rbuf);
        resize(n[j-1], a[j-1],n[j], a[j],rbuf);
        
        param[j][0] = param0[0]*(double)n[j][0]/(double)n0[0];
        param[j][1] = param0[1]*(double)n[j][1]/(double)n0[1];
        param[j][2] = param0[2]*(double)n[j][2]/(double)n0[2];
        param[j][3] = param0[3];
        param[j][4] = param0[4];
    }

    u[ng-1][0] = bo[ng-1][0]/(a[ng-1][0] + param0[4]);

    for(j=ng-2; j>=0; j--)
    {
        int jc;
        resize(n[j+1],u[j+1],n[j],u[j],rbuf);
        if(j>0) copy(m[j],bo[j],b[j]);
        for(jc=0; jc<c; jc++)
        {
            int jj;
            for(jj=j; jj<ng-1; jj++)
            {
                relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
                Atimesp(n[jj], a[jj], param[jj], u[jj], res);
                for(i=0; i<m[jj]; i++)
                    res[i] = b[jj][i] - res[i];
                resize(n[jj],res,n[jj+1],b[jj+1],rbuf);
                zeros(m[jj+1],u[jj+1]);
            }
            u[ng-1][0] = b[ng-1][0]/(a[ng-1][0] + param0[4]);
            for(jj=ng-2; jj>=j; jj--)
            {
                resize(n[jj+1],u[jj+1],n[jj],res,rbuf);
                addto(m[jj], u[jj], res);
                relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
            }
        }
    }
}

