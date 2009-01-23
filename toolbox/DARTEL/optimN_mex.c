/* $Id: optimN_mex.c 2644 2009-01-23 13:01:50Z john $ */
/* (c) John Ashburner (2007) */

#include "mex.h"
#include <math.h>
#include "optimN.h"

static void fmg_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd, i;
     int dm[4];
    int   cyc=1, nit=1, rtype=0;
    float *A, *b, *x, *scratch;
    static double param[6] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    double scal[256];

    if ((nrhs!=3 && nrhs!=4) || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and single");
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsSingle(prhs[1]))
        mexErrMsgTxt("Data must be numeric, real, full and single");

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    for(i=0; i<nd; i++) dm[i] = mxGetDimensions(prhs[1])[i];
    for(i=nd; i<4; i++) dm[i] = 1;

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    if ((nd==4) && (mxGetDimensions(prhs[0])[3] != (dm[3]*(dm[3]+1))/2))
        mexErrMsgTxt("Incompatible 4th dimension (must be (n*(n+1))/2).");
    if (nd>3) nd=3;
    for(i=0; i<nd; i++) if (mxGetDimensions(prhs[0])[i] != dm[i]) mexErrMsgTxt("Incompatible dimensions.");
    for(i=nd; i<3; i++) if (dm[i] != 1) mexErrMsgTxt("Incompatible dimensions.");
 
    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Data must be numeric, real, full and double");

    if (mxGetNumberOfElements(prhs[2]) != 9)
        mexErrMsgTxt("Third argument should contain rtype, vox1, vox2, vox3, param1, param2, param3, ncycles and relax-its.");
    rtype    = (int)(mxGetPr(prhs[2])[0]);
    param[0] = 1/mxGetPr(prhs[2])[1];
    param[1] = 1/mxGetPr(prhs[2])[2];
    param[2] = 1/mxGetPr(prhs[2])[3];
    param[3] = mxGetPr(prhs[2])[4];
    param[4] = mxGetPr(prhs[2])[5];
    param[5] = mxGetPr(prhs[2])[6];
    cyc      = mxGetPr(prhs[2])[7];
    nit      = (int)(mxGetPr(prhs[2])[8]);

    if (nrhs==4)
    {
        double *s;
        if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[3]) != dm[3])
            mexErrMsgTxt("Incompatible number of scales.");
        s = (double *)mxGetPr(prhs[3]);
        for(i=0; i< dm[3]; i++)
            scal[i] = s[i];
    }
    else
    {
        for(i=0; i<dm[3]; i++)
            scal[i] = 1.0;
    }
    plhs[0] = mxCreateNumericArray(4,(unsigned int *)dm, mxSINGLE_CLASS, mxREAL);

    A       = (float *)mxGetPr(prhs[0]);
    b       = (float *)mxGetPr(prhs[1]);
    x       = (float *)mxGetPr(plhs[0]);
    scratch = (float *)mxCalloc(fmg_scratchsize((int *)dm),sizeof(float));
    fmg((int *)dm, A, b, rtype, param, scal, cyc, nit, x, scratch);
    mxFree((void *)scratch);
}

static void vel2mom_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd, i;
     int dm[4];
    int rtype = 0;
    static double param[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    double scal[256];

    if ((nrhs!=2 && nrhs!=3) || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and single");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    for(i=0; i<nd; i++) dm[i] = mxGetDimensions(prhs[0])[i];
    for(i=nd; i<4; i++) dm[i] = 1;

    if (mxGetNumberOfElements(prhs[1]) != 7)
        mexErrMsgTxt("Parameters should contain rtype, vox1, vox2, vox3, param1, param2 and param3.");
    rtype    = (int)(mxGetPr(prhs[1])[0]);
    param[0] = 1/mxGetPr(prhs[1])[1];
    param[1] = 1/mxGetPr(prhs[1])[2];
    param[2] = 1/mxGetPr(prhs[1])[3];
    param[3] = mxGetPr(prhs[1])[4];
    param[4] = mxGetPr(prhs[1])[5];
    param[5] = mxGetPr(prhs[1])[6];

    if (nrhs==3)
    {
        double *s;
        if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[2]) != dm[3])
            mexErrMsgTxt("Incompatible number of scales.");
        s = (double *)mxGetPr(prhs[2]);
        for(i=0; i< dm[3]; i++)
            scal[i] = s[i];
    }
    else
    {
        for(i=0; i<dm[3]; i++)
            scal[i] = 1.0;
    }
 
    plhs[0] = mxCreateNumericArray(nd, (unsigned int *)dm, mxSINGLE_CLASS, mxREAL);

    if (rtype==1)
        LtLf_me((int *)dm, (float *)mxGetPr(prhs[0]), param, scal, (float *)mxGetPr(plhs[0]));
    else if (rtype==2)
        LtLf_be((int *)dm, (float *)mxGetPr(prhs[0]), param, scal, (float *)mxGetPr(plhs[0]));
    else
        mexErrMsgTxt("Regularisation type not recognised.");
}

#include<string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ((nrhs>=1) && mxIsChar(prhs[0]))
    {
        int buflen;
        char *fnc_str;
        buflen = mxGetNumberOfElements(prhs[0]);
        fnc_str = (char *)mxCalloc(buflen+1,sizeof(mxChar));
        mxGetString(prhs[0],fnc_str,buflen+1);
        if (!strcmp(fnc_str,"vel2mom"))
        {
            mxFree(fnc_str);
            vel2mom_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"fmg")  || !strcmp(fnc_str,"FMG"))
        {
            mxFree(fnc_str);
            fmg_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else
        {
            mxFree(fnc_str);
            mexErrMsgTxt("Option not recognised.");
        }
    }
    else
    {
        fmg_mexFunction(nlhs, plhs, nrhs, prhs);
    }
}

