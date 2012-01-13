/*
 * $Id: spm_bsplins.c 4624 2012-01-13 13:27:08Z john $
 * John Ashburner
 */

/*
 * This code is based on that of Philippe Thevenaz, which I took from:
 *  http://bigwww.epfl.ch/algorithms.html
 *
 * It has been substantially modified, so blame me (John Ashburner) if there
 * are any bugs. Many thanks to Philippe Thevenaz for advice with the code.
 *
 * See:
 *  M. Unser.
 *  "Splines: A Perfect Fit for Signal and Image Processing,"
 *  IEEE Signal Processing Magazine, 16(6):22-38 (1999)
 *
 *  P. Thevenaz and T. Blu and M. Unser.
 *  "Interpolation Revisited"
 *  IEEE Transactions on Medical Imaging 19(7):739-758 (2000).
*/

#include <math.h>
#include "mex.h"
#include "bsplines.h"

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

static void fun(double c[], int m0, int m1, int m2,
    int n, double x0[], double x1[], double x2[], int d[],
    int cond, int (*bnd[])(), double f[])
{
    int j;
    double NaN = mxGetNaN();

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
static void dfun(double c[], int m0, int m1, int m2,
    int n, double x0[], double x1[], double x2[],int d[],
    int cond, int (*bnd[])(),
    double f[], double df0[], double df1[], double df2[])
{
    int j;
    double NaN = mxGetNaN();

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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int k, d[3], n, nd;
    int m0=1, m1=1, m2=1;
    double *x0, *x1, *x2, *c, *f, *df0, *df1, *df2;
    const mwSize *dims;
    int (*bnd[3])();
    int cond;

    /* Usage:
            f = function(c,x0,x1,x2,d)
                c - B-spline coefficients
                x0, x1, x2 - co-ordinates
                d   - B-spline degree
                f   - sampled function
       or:
            [f,df0,df1,df2] = function(c,x0,x1,x2,d)
                c - B-spline coefficients
                x0, x1, x2 - co-ordinates
                d   - B-spline degree
                f   - sampled function
                df0, df1, df2   - sampled derivatives
    */
    if (nrhs < 5 || nlhs>4)
        mexErrMsgTxt("Incorrect usage.");

    for(k=0; k<5; k++)
    {
        if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
            mxIsSparse(prhs[k]) || !mxIsDouble(prhs[k]))
            mexErrMsgTxt("Input must be numeric, real, full and double precision.");
    }

    if ((mxGetM(prhs[4])*mxGetN(prhs[4]) != 3) && (mxGetM(prhs[4])*mxGetN(prhs[4]) != 6))
        mexErrMsgTxt("Incorrect usage.");

    /* Degree of spline */
    for(k=0; k<3; k++)
    {
        d[k] = floor(mxGetPr(prhs[4])[k]+0.5);
        if (d[k]<0 || d[k]>7)
            mexErrMsgTxt("Bad spline degree.");
    }

    cond = 0;
    for(k=0; k<3; k++) bnd[k] = mirror;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) == 6)
    {
        for(k=0; k<3; k++)
            if (mxGetPr(prhs[4])[k+3])
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
    nd = mxGetNumberOfDimensions(prhs[1]);
    dims = mxGetDimensions(prhs[1]);
    if (mxGetNumberOfDimensions(prhs[2]) != nd || mxGetNumberOfDimensions(prhs[3]) != nd)
        mexErrMsgTxt("Incompatible dimensions.");
    n = 1;
    for(k=0; k<nd; k++)
    {
        if (mxGetDimensions(prhs[2])[k] != dims[k] || mxGetDimensions(prhs[3])[k] != dims[k])
            mexErrMsgTxt("Incompatible dimensions.");
        n *=dims[k];
    }

    /* Sampled data same size as sampling co-ords */
    plhs[0] = mxCreateNumericArray(nd,dims, mxDOUBLE_CLASS, mxREAL);

    /* Pointers to double precision data */
    c  = mxGetPr(prhs[0]);
    x0 = mxGetPr(prhs[1]);
    x1 = mxGetPr(prhs[2]);
    x2 = mxGetPr(prhs[3]);
    f  = mxGetPr(plhs[0]);

    if (nlhs<=1)
        fun(c, m0,m1,m2, n, x0,x1,x2, d, cond,bnd, f);
    else
    {
        plhs[1] = mxCreateNumericArray(nd,dims, mxDOUBLE_CLASS, mxREAL);
        plhs[2] = mxCreateNumericArray(nd,dims, mxDOUBLE_CLASS, mxREAL);
        plhs[3] = mxCreateNumericArray(nd,dims, mxDOUBLE_CLASS, mxREAL);
        df0 = mxGetPr(plhs[1]);
        df1 = mxGetPr(plhs[2]);
        df2 = mxGetPr(plhs[3]);
        dfun(c, m0,m1,m2, n, x0,x1,x2, d, cond,bnd, f,df0,df1,df2);
    }
}
