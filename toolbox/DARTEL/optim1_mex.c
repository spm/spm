/* $Id: optim1_mex.c 1137 2008-02-06 15:58:21Z spm $ */
/* (c) John Ashburner (2007) */

#include "mex.h"
#include <math.h>
#include "optim1.h"

void cgs_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int *dm;
    int        nit=1, rtype=2;
    double      tol=1e-10, *A, *b, *x, *scratch1, *scratch2, *scratch3;
    static double param[] = {1.0, 1.0, 1.0, 1.0, 0.0};

    if (nrhs!=3 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfDimensions(prhs[0])!=3) mexErrMsgTxt("Wrong number of dimensions.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfDimensions(prhs[1])!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[1]);

    if (mxGetDimensions(prhs[0])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[0])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");
    if (mxGetDimensions(prhs[0])[2] != dm[2])
        mexErrMsgTxt("Incompatible 3rd dimension.");

    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfElements(prhs[2]) != 4)
        mexErrMsgTxt("Third argument should contain rtype, lambda, tol and nit.");

    rtype    = (int)(mxGetPr(prhs[2])[0]);
    param[3] = mxGetPr(prhs[2])[1];
    tol      = mxGetPr(prhs[2])[2];
    nit      = (int)(mxGetPr(prhs[2])[3]);

    plhs[0] = mxCreateNumericArray(3,dm, mxDOUBLE_CLASS, mxREAL);

    A       = (double *)mxGetPr(prhs[0]);
    b       = (double *)mxGetPr(prhs[1]);
    x       = (double *)mxGetPr(plhs[0]);

    scratch1 = (double *)mxCalloc(dm[0]*dm[1]*dm[2],sizeof(double));
    scratch2 = (double *)mxCalloc(dm[0]*dm[1]*dm[2],sizeof(double));
    scratch3 = (double *)mxCalloc(dm[0]*dm[1]*dm[2],sizeof(double));

    cgs((int *)dm, A, b, rtype, param, tol, nit, x,scratch1,scratch2,scratch3);

    mxFree((void *)scratch3);
    mxFree((void *)scratch2);
    mxFree((void *)scratch1);
}

void fmg_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int *dm;
    int        cyc=1, nit=1, rtype=2;
    double     *A, *b, *x, *scratch;
    static double param[] = {1.0, 1.0, 1.0, 1.0, 0.0};

    if (nrhs!=3 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfDimensions(prhs[0])!=3) mexErrMsgTxt("Wrong number of dimensions.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfDimensions(prhs[1])!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[1]);

    if (mxGetDimensions(prhs[0])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[0])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");
    if (mxGetDimensions(prhs[0])[2] != dm[2])
        mexErrMsgTxt("Incompatible 3rd dimension.");

    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfElements(prhs[2]) != 4)
        mexErrMsgTxt("Third argument should contain rtype, lambda, ncycles and relax-its.");

    rtype    = (int)(mxGetPr(prhs[2])[0]);
    param[3] = mxGetPr(prhs[2])[1];
    cyc      = (int)(mxGetPr(prhs[2])[2]);
    nit      = (int)(mxGetPr(prhs[2])[3]);

    plhs[0] = mxCreateNumericArray(3,dm, mxDOUBLE_CLASS, mxREAL);

    A       = (double *)mxGetPr(prhs[0]);
    b       = (double *)mxGetPr(prhs[1]);
    x       = (double *)mxGetPr(plhs[0]);
    scratch = (double *)mxCalloc(fmg_scratchsize((int *)dm),sizeof(double));
    fmg((int *)dm, A, b, rtype, param, cyc, nit, x, scratch);
    mxFree((void *)scratch);
}


void rsz_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nc[3];
    const int *na;
    double *a, *b, *c;
    if ((nrhs!=2) || (nlhs>1))
        mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
            mexErrMsgTxt("Data must be numeric, real, full and double");

    if (mxGetNumberOfDimensions(prhs[0])!=3) mexErrMsgTxt("Wrong number of dimensions.");
    na = mxGetDimensions(prhs[0]);

    if (mxGetNumberOfElements(prhs[1]) != 3)
    {
        mexErrMsgTxt("Dimensions argument is wrong size.");
    }
    nc[0] = (int)mxGetPr(prhs[1])[0];
    nc[1] = (int)mxGetPr(prhs[1])[1];
    nc[2] = (int)mxGetPr(prhs[1])[2];

    a = (double *)mxGetPr(prhs[0]);
    b = (double *)mxCalloc(na[0]*nc[1]+3*nc[0]*nc[1],sizeof(double));
    plhs[0] = mxCreateNumericArray(3,nc, mxDOUBLE_CLASS, mxREAL);
    c = (double *)mxGetPr(plhs[0]);
    resize((int *)na, a, nc, c, b);
    (void)mxFree(b);
}

void LLvbe_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd;
    const int *dm;
    static double param[] = {1.0, 1.0, 1.0, 1.0, 0.0};

    if (nrhs!=1 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);

    plhs[0] = mxCreateNumericArray(nd,dm, mxDOUBLE_CLASS, mxREAL);
    LtLf_be((int *)dm, (double *)mxGetPr(prhs[0]), param, (double *)mxGetPr(plhs[0]));
}

void LLvme_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd;
    const int *dm;
    static double param[] = {1.0, 1.0, 1.0, 1.0, 0.0};

    if (nrhs!=1 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);

    plhs[0] = mxCreateNumericArray(nd,dm, mxDOUBLE_CLASS, mxREAL);
    LtLf_me((int *)dm, (double *)mxGetPr(prhs[0]), param, (double *)mxGetPr(plhs[0]));
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
        if (!strcmp(fnc_str,"LLvbe"))
        {
            mxFree(fnc_str);
            LLvbe_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"LLvme"))
        {
            mxFree(fnc_str);
            LLvme_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"fmg")  || !strcmp(fnc_str,"FMG"))
        {
            mxFree(fnc_str);
            fmg_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"cgs")  || !strcmp(fnc_str,"CGS"))
        {
            mxFree(fnc_str);
            cgs_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"rsz")  || !strcmp(fnc_str,"RSZ"))
        {
            mxFree(fnc_str);
            rsz_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
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

