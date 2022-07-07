/*
 * John Ashburner
 * Copyright (C) 1999-2022 Wellcome Centre for Human Neuroimaging
 */

#include <math.h>
#include "mex.h"
#include "hist2.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mwSize *dimsp;
    mwSize dimsf[3], dimsg[3];
    mwSize nd, i;
    float s[3];

    if (nrhs>4 || nrhs<3 || nlhs>1) mexErrMsgTxt("Incorrect usage.");

    if (!mxIsNumeric(prhs[0]) || !mxIsUint8(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("Wrong sort of data (1).");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>3) mexErrMsgTxt("Wrong number of dims (1).");
    dimsp = mxGetDimensions(prhs[0]);
    for(i=0; i<nd; i++) dimsg[i] = dimsp[i];
    for(i=nd; i<3; i++) dimsg[i] = 1;


    if (!mxIsNumeric(prhs[1]) || !mxIsUint8(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgTxt("Wrong sort of data (2).");
    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd>3) mexErrMsgTxt("Wrong number of dims (2).");
    dimsp = mxGetDimensions(prhs[1]);
    for(i=0; i<nd; i++) dimsf[i] = dimsp[i];
    for(i=nd; i<3; i++) dimsf[i] = 1;


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
