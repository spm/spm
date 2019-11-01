/*
 * $Id: shoot_bsplines.c 7685 2019-11-01 12:56:19Z john $
 * John Ashburner
 */
 
/*
 * This code is a modified version of that of Philippe Thevenaz, which I took from:
 *  http://bigwww.epfl.ch/algorithms.html
 *
 * It has been substantially modified, so blame me (John Ashburner) if there
 * are any bugs. Many thanks to Philippe Thevenaz for advice with the code.
 *
 * See:
 *  M. Unser, A. Aldroubi and M. Eden.
 *  "B-Spline Signal Processing: Part I-Theory,"
 *  IEEE Transactions on Signal Processing 41(2):821-832 (1993).
 *
 *  M. Unser, A. Aldroubi and M. Eden.
 *  "B-Spline Signal Processing: Part II-Efficient Design and Applications,"
 *  IEEE Transactions on Signal Processing 41(2):834-848 (1993).
 *
 *  M. Unser.
 *  "Splines: A Perfect Fit for Signal and Image Processing,"
 *  IEEE Signal Processing Magazine 16(6):22-38 (1999).
 *
*/

#include <math.h>
#include "mex.h"
#include "bsplines.h"


/***************************************************************************************
Deconvolve the B-spline basis functions from the image volume
    vol - a handle for the volume to deconvolve
    c - the coefficients (arising from the deconvolution)
    d - the spline degree
    splinc0, splinc1, splinc2   - functions for 1D deconvolutions
*/
static int vol_coeffs(mwSize vdim[], float vol[], float c[], int d[], void (*splinc[])())
{
    double  p[4];
    float *cp;
    int np;
    mwSize i, j, k, n;
    float f[10240];

    /* Check that dimensions don't exceed size of f */
    if (vdim[1]>10240 ||vdim[2]>10240)
        return(1);

    /* Do a straight copy */
    cp = c;
    for(cp=c; cp<c+vdim[2]*vdim[1]*vdim[0]; cp++, vol++)
    {
        *cp = *vol;
        if (!mxIsFinite(*cp)) *cp = 0.0;

    }

    /* Deconvolve along the fastest dimension (X) */
    if (d[0]>1 && vdim[0]>1)
    {
        if (get_poles(d[0], &np, p)!=0) return(1);
#       pragma omp parallel for collapse(2) private(cp)
        for(k=0; k<vdim[2]; k++)
        {
            /* double dk = k+1; */
            for(j=0; j<vdim[1]; j++)
            {
                cp = &c[vdim[0]*(j+vdim[1]*k)];
                splinc[0](cp, vdim[0], p, np);
            }
        }
    }

    /* Deconvolve along the middle dimension (Y) */
    if (d[1]>1 && vdim[1]>1)
    {
        if (get_poles(d[1], &np, p)!=0) return(1);
        n =vdim[0];
#       pragma omp parallel for collapse(2) private(cp,f)
        for(k=0; k<vdim[2]; k++)
        {
            for(i=0;i<vdim[0];i++)
            {
                cp = &c[i+vdim[0]*vdim[1]*k];
                for(j=0; j<vdim[1]; j++, cp+=n) f[j] = *cp;
                splinc[1](f, vdim[1], p, np);
                cp = &c[i+vdim[0]*vdim[1]*k];
                for(j=0; j<vdim[1]; j++, cp+=n) *cp = f[j];
            }
        }
    }

    /* Deconvolve along the slowest dimension (Z) */
    if (d[2]>1 && vdim[2]>1)
    {
        if (get_poles(d[2], &np, p)!=0) return(1);
        n = vdim[0]*vdim[1];
#       pragma omp parallel for collapse(2) private(cp,f)
        for(j=0; j<vdim[1]; j++)
        {
            for(i=0;i<vdim[0];i++)
            {
                cp = &c[i+vdim[0]*j];
                for(k=0; k<vdim[2]; k++, cp+=n) f[k] = *cp;
                splinc[2](f, vdim[2], p, np);
                cp = &c[i+vdim[0]*j];
                for(k=0; k<vdim[2]; k++, cp+=n) *cp = f[k];
            }
        }
    }
    return(0);
}

/***************************************************************************************
*/
void bsplinc_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize vdim[3];
    int k, nd, d[3];
    float *c, *f;
    void (*splinc[3])();

    if (nrhs < 2 || nlhs > 1)
        mexErrMsgTxt("Incorrect usage.");
    if (mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) ||
        (mxGetM(prhs[1])*mxGetN(prhs[1]) != 3 && mxGetM(prhs[1])*mxGetN(prhs[1]) != 6))
        mexErrMsgTxt("Incorrect usage.");

    for(k=0; k<3; k++)
    {
        d[k] = floor(mxGetPr(prhs[1])[k]+0.5);
        if (d[k]<0 || d[k]>7)
            mexErrMsgTxt("Bad spline degree.");
    }

    for(k=0; k<3; k++) splinc[k] = splinc_mirror;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) == 6)
    {
        for(k=0; k<3; k++)
            if (mxGetPr(prhs[1])[k+3])
                splinc[k] = splinc_wrap;
    }

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ||
         mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
         mexErrMsgTxt("Input must be numeric, real, full and single precision.");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>3) mexErrMsgTxt("Too many dimensions.");
    for(k=0; k<nd; k++) vdim[k] =  mxGetDimensions(prhs[0])[k];
    for(k=nd; k<3; k++) vdim[k] = 1;
    f = (float *)mxGetPr(prhs[0]);

    plhs[0] = mxCreateNumericArray(3,vdim, mxSINGLE_CLASS, mxREAL);
    c = (float *)mxGetPr(plhs[0]);

    vol_coeffs(vdim, f, c, d, splinc);
}

/***************************************************************************************
Loop through data and resample the points
    c   - Volume of B-spline coefficients
    m0,m1,m2    - dimensions of c
    n   - number of points to resample
    x0,x1,x2    - array of co-ordinate to sample
    d   - degree of spline used
    cond    - code determining boundaries to mask at
    bnd - functions for dealing with edges
    f   - resampled data
*/
#define TINY 5e-2

static void fun(float c[], int m0, int m1, int m2,
    mwSize n, float x0[], float x1[], float x2[], int d[],
    int cond, int (*bnd[])(), float f[])
{
    mwSize j;
    double NaN = mxGetNaN();

#   pragma omp parallel for
    for(j=0; j<n; j++)
    {
        if (((cond&1) | (x0[j]>=1-TINY && x0[j]<=m0+TINY)) &&
            ((cond&2) | (x1[j]>=1-TINY && x1[j]<=m1+TINY)) &&
            ((cond&4) | (x2[j]>=1-TINY && x2[j]<=m2+TINY)))
            f[j] = sample3(c, m0,m1,m2, x0[j]-1,x1[j]-1,x2[j]-1, d, bnd);
        else
            f[j] = NaN;
    }
}


/***************************************************************************************
Loop through data and resample the points and their derivatives
    c   - Volume of B-spline coefficients
    m0,m1,m2    - dimensions of c
    n   - number of points to resample
    x0,x1,x2    - array of co-ordinate to sample
    d   - degrees of splines used
    cond    - code determining boundaries to mask at
    bnd - functions for dealing with edges
    f   - resampled data
    df0, df1, df2   - gradients
*/
static void dfun(float c[], int m0, int m1, int m2,
    mwSize n, float x0[], float x1[], float x2[], int d[],
    int cond, int (*bnd[])(),
    float f[], float df0[], float df1[], float df2[])
{
    mwSize j;
    float NaN = mxGetNaN();

    for(j=0; j<n; j++)
    {
        if (((cond&1) | (x0[j]>=1-TINY && x0[j]<=m0+TINY)) &&
            ((cond&2) | (x1[j]>=1-TINY && x1[j]<=m1+TINY)) &&
            ((cond&4) | (x2[j]>=1-TINY && x2[j]<=m2+TINY)))
            f[j] = dsample3(c, m0,m1,m2, x0[j]-1,x1[j]-1,x2[j]-1, d,
                &df0[j],&df1[j],&df2[j], bnd);
        else
            f[j] = NaN;
    }
}


/***************************************************************************************/
void bsplins_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int k, nd, d[3];
    mwSize n;
    int m0=1, m1=1, m2=1;
    float *x0, *x1, *x2, *c, *f, *df0, *df1, *df2;
    const mwSize *dims;
    int (*bnd[3])();
    int cond;

    /* Usage:
            f = function(c,Y,d)
                c - B-spline coefficients
                Y - co-ordinates
                d - B-spline degree
                f - sampled function
       or:
            [f,df0,df1,df2] = function(c,Y,d)
                c - B-spline coefficients
                Y - co-ordinates
                d - B-spline degree
                f - sampled function
                df0, df1, df2   - sampled derivatives
    */
    if (nrhs < 3 || nlhs>4)
        mexErrMsgTxt("Incorrect usage.");

    for(k=0; k<2; k++)
    {
        if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
            mxIsSparse(prhs[k]) || !mxIsSingle(prhs[k]))
            mexErrMsgTxt("Input must be numeric, real, full and single precision.");
    }

    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) ||
        mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Input must be numeric, real, full and double precision.");

    if ((mxGetM(prhs[2])*mxGetN(prhs[2]) != 3) && (mxGetM(prhs[2])*mxGetN(prhs[2]) != 6))
        mexErrMsgTxt("Incorrect usage.");

    /* Degree of spline */
    for(k=0; k<3; k++)
    {
        d[k] = floor(mxGetPr(prhs[2])[k]+0.5);
        if (d[k]<0 || d[k]>7)
            mexErrMsgTxt("Bad spline degree.");
    }

    cond = 0;
    for(k=0; k<3; k++) bnd[k] = mirror;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) == 6)
    {
        for(k=0; k<3; k++)
            if (mxGetPr(prhs[2])[k+3])
            {
                bnd[k] = wrap;
                cond += 1<<k;
            }
    }

    /* if (d==0 && nlhs>1)
        mexErrMsgTxt("Cant compute gradients when using B-spline(0) interp."); */

    /* Dimensions of coefficient volume */
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>3) mexErrMsgTxt("Too many coefficient dimensions.");
    dims = mxGetDimensions(prhs[0]);
    if (nd>=1) m0 = dims[0];
    if (nd>=2) m1 = dims[1];
    if (nd>=3) m2 = dims[2];

    /* Dimensions of sampling co-ordinates */
    nd   = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dims = mxGetDimensions(prhs[1]);
    n = 1;
    for(k=0; k<3; k++) n *=dims[k];

    /* Sampled data same size as sampling co-ords */
    plhs[0] = mxCreateNumericArray(3,dims, mxSINGLE_CLASS, mxREAL);

    /* Pointers to double precision data */
    c  = (float *)mxGetPr(prhs[0]);
    x0 = (float *)mxGetPr(prhs[1]);
    x1 = (float *)mxGetPr(prhs[1])+n;
    x2 = (float *)mxGetPr(prhs[1])+n*2;
    f  = (float *)mxGetPr(plhs[0]);

    if (nlhs<=1)
        fun(c, m0,m1,m2, n, x0,x1,x2, d, cond,bnd, f);
    else
    {
        plhs[1] = mxCreateNumericArray(3,dims, mxSINGLE_CLASS, mxREAL); df0 = (float *)mxGetPr(plhs[1]);
        plhs[2] = mxCreateNumericArray(3,dims, mxSINGLE_CLASS, mxREAL); df1 = (float *)mxGetPr(plhs[2]);
        plhs[3] = mxCreateNumericArray(3,dims, mxSINGLE_CLASS, mxREAL); df2 = (float *)mxGetPr(plhs[3]);
        dfun(c, m0,m1,m2, n, x0,x1,x2, d, cond,bnd, f,df0,df1,df2);
    }
}

