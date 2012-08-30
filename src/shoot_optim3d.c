/* $Id: shoot_optim3d.c 4875 2012-08-30 20:04:30Z john $ */
/* (c) John Ashburner (2011) */

#include<mex.h>
#include<math.h>
extern double log(double x);

#include "shoot_optim3d.h"
#include "shoot_multiscale.h"
#include "shoot_regularisers.h"

static double dotprod(mwSize m, float a[], float b[])
{
    mwSize i;
    double dp = 0.0;
    for(i=0; i<m; i++)
        dp += a[i]*b[i];
    return(dp);
}

static void addscaled(mwSize m, float a[], float b[], double s)
{
    mwSize i;
    for(i=0; i<m; i++)
        a[i] += s*b[i];
}

float norm(mwSize m, float a[])
{
    mwSize i;
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

void cgs3(mwSize dm[], float A[], float b[], double param[], double tol, int nit,
             float x[], float r[], float p[], float Ap[])
{
    mwSize i, m = dm[0]*dm[1]*dm[2]*3, it;
    double rtr, nb, rtrold, alpha, beta;

    /* printf("\n **** %dx%d ****\n",dm[0],dm[1]); */

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

static void rescale(mwSize n, float *a, double s)
{
    mwSize i;
    for(i=0; i<n; i++)
        a[i] *= s;
}

static void restrict_g(mwSize na[], float *a, float *c, float *b)
{
    mwSize i, nc[3], m;
    for(i=0; i<3; i++)
        nc[i] = ceil(na[i]/2.0);
    m = nc[0]*nc[1]*nc[2];
    for(i=0; i<3; i++)
    {
        restrict_vol(na, a+i*na[0]*na[1]*na[2], nc, c+i*m, b);
        rescale(m,c+i*m,(double)na[i]/(double)nc[i]);
    }
}

static void restrict_h(mwSize na[], float *a, float *c, float *b)
{
    mwSize i, nc[3], m;
    double s[3];
    for(i=0; i<3; i++)
    {
        nc[i] = ceil(na[i]/2.0);
        s[i]  = (double)na[i]/(double)nc[i];
    }
    m = nc[0]*nc[1]*nc[2];

    for(i=0; i<6; i++)
        restrict_vol(na, a+i*na[0]*na[1]*na[2], nc, c+i*m, b);

    for(i=0; i<3; i++)
        rescale(m,c+i*m,s[i]*s[i]);
    rescale(m,c+3*m,s[0]*s[1]);
    rescale(m,c+4*m,s[0]*s[2]);
    rescale(m,c+5*m,s[1]*s[2]);
}

static void prolong(mwSize na[], float *a, mwSize nc[], float *c, float *b)
{
    mwSize i,m = nc[0]*nc[1]*nc[2];
    for(i=0; i<3; i++)
    {
        resize_vol(na, a+i*na[0]*na[1]*na[2], nc, c+i*m, b);
        rescale(m,c+i*m,(double)nc[i]/(double)na[i]);
    }
}

static void zeros(mwSize n, float *a)
{
    mwSize i;
    for(i=0; i<n; i++)
        a[i] = 0.0;
}

static void copy(mwSize n, float *a, float *b)
{
    mwSize i;
    for(i=0; i<n; i++)
        b[i] = a[i];
}

static void addto(mwSize n, float *a, float *b)
{
    mwSize i;
    for(i=0; i<n; i++)
        a[i] += b[i];
}

mwSize fmg3_scratchsize(mwSize n0[], int use_hessian)
{
    mwSize n[64][3], m[64], bs, j, num_blocks;
    if (use_hessian)
        num_blocks = 15; /* Uses a further 6 volumes to represent a symmetric tensor field */
    else
        num_blocks = 9;

    /* Figure out bs, which is the sum of the number of elements
       of a 3D volume over all scales (except the highest). */
    bs = 0;
    n[0][0] = n0[0];
    n[0][1] = n0[1];
    n[0][2] = n0[2];
    for(j=1; j<64; j++)
    {
        n[j][0] = ceil(n[j-1][0]/2.0);
        n[j][1] = ceil(n[j-1][1]/2.0);
        n[j][2] = ceil(n[j-1][2]/2.0);
        m[j]    = n[j][0]*n[j][1]*n[j][2];
        bs += m[j];
        if ((n[j][0]<2) && (n[j][1]<2) && (n[j][2]<2))
            break;
    }
    return((3*n0[0]*n0[1]*n0[2] + n[0][0]*n[1][1]+4*n[0][0]*n[0][1] + num_blocks*bs));
}

/*
    Full Multigrid solver.  See Numerical Recipes (second edition) for more
    information
*/
void fmg3(mwSize n0[], float *a0, float *b0, double param0[], int c, int nit,
          float *u0, float *scratch)
{
    mwSignedIndex i, j, ng, bs;
    mwSize n[64][3], m[64];
    float *bo[64], *a[64], *b[64], *u[64], *res, *rbuf;
    double param[64][8];

    /* Dimensions of native resolution grids */
    n[0][0] = n0[0];
    n[0][1] = n0[1];
    n[0][2] = n0[2];
    m[0]    = n0[0]*n0[1]*n0[2];

    /* Dimensions of lower resolution grids */
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

    /* Set up pointers to native data */
    bo[0]   = b0;
    b[0]    = b0;
    u[0]    = u0;
    a[0]    = a0;

    /* Set up pointers to scratch space */
    res          = scratch;           /* Residuals (defect) */
    rbuf         = scratch + 3*m[0];  /* Additional memory needed for restriction/prolongation */

    bo[1]        = scratch + 3*m[0] + n[0][0]*n[1][1]+4*n[0][0]*n[0][1];
    b[1]         = scratch + 3*m[0] + n[0][0]*n[1][1]+4*n[0][0]*n[0][1] + 3*bs;
    u[1]         = scratch + 3*m[0] + n[0][0]*n[1][1]+4*n[0][0]*n[0][1] + 6*bs;
    if (a0) a[1] = scratch + 3*m[0] + n[0][0]*n[1][1]+4*n[0][0]*n[0][1] + 9*bs;
    else    a[1] = 0;
    for(j=2; j<ng; j++)
    {
        bo[j]    = bo[j-1]+3*m[j-1];
        b[j]     =  b[j-1]+3*m[j-1];
        u[j]     =  u[j-1]+3*m[j-1];
        if (a0)
            a[j] =  a[j-1]+6*m[j-1];
        else
            a[j] = 0;
    }

    /* Create grids of gradients and hessians, as well as adjusting
       the (reciprocals of the) voxel sizes used for defining the
       regularisation operators */
    param[0][0] = param0[0]; /* 1/vox_x */
    param[0][1] = param0[1]; /* 1/vox_y */
    param[0][2] = param0[2]; /* 1/vox_z */
    param[0][3] = param0[3]; /* 1st regularisation parameter */
    param[0][4] = param0[4]; /* 2nd regularisation parameter */
    param[0][5] = param0[5]; /* 3rd regularisation parameter */
    param[0][6] = param0[6]; /* 4th regularisation parameter */
    param[0][7] = param0[7]; /* 5th regularisation parameter */

    for(j=1; j<ng; j++)
    {
        restrict_g(n[j-1],bo[j-1],bo[j],rbuf);
        if (a0)
            restrict_h(n[j-1],a[j-1],a[j],rbuf);

        param[j][0] = param0[0]*(double)n[j][0]/n0[0];
        param[j][1] = param0[1]*(double)n[j][1]/n0[1];
        param[j][2] = param0[2]*(double)n[j][2]/n0[2];
        param[j][3] = param[0][3];
        param[j][4] = param[0][4];
        param[j][5] = param[0][5];
        param[j][6] = param[0][6];
        param[j][7] = param[0][7];
    }

    if (u[0][0]==0) /* No starting estimate so do Full Multigrid */
    {
        relax(n[ng-1], a[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]); 
        for(j=ng-2; j>=0; j--)
        {
            mwSignedIndex jc;
            prolong(n[j+1],u[j+1],n[j],u[j],rbuf);
            if(j>0) copy(3*m[j],bo[j],b[j]);
            for(jc=0; jc<c; jc++)
            {
                mwSignedIndex jj;
                for(jj=j; jj<ng-1; jj++) /* From high res to lowest res */
                {
                    relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
                    Atimesp(n[jj], a[jj], param[jj], u[jj], res);
                    for(i=0; i<3*m[jj]; i++)
                        res[i] = b[jj][i] - res[i];
                    restrict_g(n[jj],res,b[jj+1],rbuf);
                    zeros(3*m[jj+1],u[jj+1]);
                }
                relax(n[ng-1], a[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]); 
                for(jj=ng-2; jj>=j; jj--) /* From lowest res to high res */
                {
                    prolong(n[jj+1],u[jj+1],n[jj],res,rbuf);
                    addto(3*m[jj], u[jj], res);
                    relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
                }
            }
        }
    }
    else /* Use starting estimate and just run some V-cycles */
    {
        mwSignedIndex jc;
/*        for(j=1; j<ng; j++)
            restrict_g(n[j-1],u[j-1],u[j],rbuf); */

        for(jc=0; jc<c; jc++) /* Loop over V-cycles */
        {
            mwSignedIndex jj;
            for(jj=0; jj<ng-1; jj++) /* From highest res to lowest res */
            {
                relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
                Atimesp(n[jj], a[jj], param[jj], u[jj], res);
                for(i=0; i<3*m[jj]; i++)
                    res[i] = b[jj][i] - res[i];
                restrict_g(n[jj],res,b[jj+1],rbuf);
                zeros(3*m[jj+1],u[jj+1]);
            }
            relax(n[ng-1], a[ng-1], b[ng-1], param[ng-1], nit, u[ng-1]); 
            for(jj=ng-2; jj>=0; jj--) /* From lowest res to highest res */
            {
                prolong(n[jj+1],u[jj+1],n[jj],res,rbuf);
                addto(3*m[jj], u[jj], res);
                relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);

            }
        }
    }
}

