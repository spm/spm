/*
 * $Id: spm_hist2.c 7629 2019-06-27 12:35:45Z john $
 * John Ashburner
 */

#include <math.h>
#include "mex.h"
#include "hist2.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mwSize *dimsf, *dimsg;
    float s[3];

    if (nrhs>4 || nrhs<3 || nlhs>1) mexErrMsgTxt("Incorrect usage.");

    if (!mxIsNumeric(prhs[0]) || !mxIsUint8(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("Wrong sort of data (1).");
    if (mxGetNumberOfDimensions(prhs[0]) != 3) mexErrMsgTxt("Wrong number of dims (1).");
    dimsg = mxGetDimensions(prhs[0]);

    if (!mxIsNumeric(prhs[1]) || !mxIsUint8(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgTxt("Wrong sort of data (2).");
    if (mxGetNumberOfDimensions(prhs[1]) != 3) mexErrMsgTxt("Wrong number of dims (2).");
    dimsf = mxGetDimensions(prhs[1]);

    if (!mxIsNumeric(prhs[2]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
        mexErrMsgTxt("Wrong sort of matrix.");
    if (mxGetM(prhs[2]) != 4 || mxGetN(prhs[2]) != 4)
        mexErrMsgTxt("Matrix must be 4x4.");

    if (nrhs == 4)
    {
        if (!mxIsNumeric(prhs[3]) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
             mxGetM(prhs[3])*mxGetN(prhs[3]) != 3)
            mexErrMsgTxt("Invalid skips.");
        s[0] = mxGetPr(prhs[3])[0];
        s[1] = mxGetPr(prhs[3])[1];
        s[2] = mxGetPr(prhs[3])[2];
    }
    else
    {
        s[0] = s[1] = s[2] = 1.0;
    }

    plhs[0] = mxCreateDoubleMatrix(256,256,mxREAL);

    hist2(mxGetPr(prhs[2]), (unsigned char *)mxGetData(prhs[0]), (unsigned char *)mxGetData(prhs[1]),
        dimsg, dimsf, mxGetPr(plhs[0]), s);

}

