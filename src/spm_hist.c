/*
 * $Id: spm_hist.c 4453 2011-09-02 10:47:25Z guillaume $
 * John Ashburner
 */

#include "mex.h"

void create_hist(mwSize n, unsigned char c[], double f[], double h[])
{
    mwIndex i;
    for(i=0; i<n; i++)
        h[c[i]] += f[i];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize n;
    if (nrhs != 2)
        mexErrMsgTxt("Incorrect usage.");

    n = mxGetM(prhs[0])*mxGetN(prhs[0]);
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != n) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsUint8(prhs[0])) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsDouble(prhs[1])) mexErrMsgTxt("Incorrect usage.");
    plhs[0] = mxCreateDoubleMatrix(256,1,mxREAL);
    create_hist(n,(unsigned char *)mxGetData(prhs[0]),mxGetPr(prhs[1]),mxGetPr(plhs[0]));
}
