/* $Id: optimizer2d_mex.c 964 2007-10-19 16:35:34Z john $ */
/* (c) John Ashburner (2007) */

#include "mex.h"
#include <math.h>
#include "optimizer2d.h"

void cgs2_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int *dm;
    int        nit=1, rtype=2;
    double     tol=1e-10, *A, *b, *x, *scratch1, *scratch2, *scratch3;
    static double param[] = {1.0, 1.0, 1.0, 0.0, 0.0};

    if (nrhs!=3 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfDimensions(prhs[0])!=3) mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetDimensions(prhs[0])[2]!=3)
        mexErrMsgTxt("3rd dimension of 1st arg must be 3.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfDimensions(prhs[1])!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[1]);
    if (dm[2]!=2)
        mexErrMsgTxt("3rd dimension of second arg must be 2.");

    if (mxGetDimensions(prhs[0])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[0])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");

    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfElements(prhs[2]) >8)
        mexErrMsgTxt("Third argument should contain rtype, vox1, vox2, param1, param2, param3, tol and nit.");
    if (mxGetNumberOfElements(prhs[2]) >=1) rtype    = (int)mxGetPr(prhs[2])[0];
    if (mxGetNumberOfElements(prhs[2]) >=2) param[0] = 1/mxGetPr(prhs[2])[1];
    if (mxGetNumberOfElements(prhs[2]) >=3) param[1] = 1/mxGetPr(prhs[2])[2];
    if (mxGetNumberOfElements(prhs[2]) >=4) param[2] = mxGetPr(prhs[2])[3];
    if (mxGetNumberOfElements(prhs[2]) >=5) param[3] = mxGetPr(prhs[2])[4];
    if (mxGetNumberOfElements(prhs[2]) >=6) param[4] = mxGetPr(prhs[2])[5];
    if (mxGetNumberOfElements(prhs[2]) >=7) tol      = mxGetPr(prhs[2])[6];
    if (mxGetNumberOfElements(prhs[2]) >=8) nit      = (int)mxGetPr(prhs[2])[7];

    plhs[0] = mxCreateNumericArray(3,dm, mxDOUBLE_CLASS, mxREAL);

    A       = mxGetPr(prhs[0]);
    b       = mxGetPr(prhs[1]);
    x       = mxGetPr(plhs[0]);

    scratch1 = (double *)mxCalloc(dm[0]*dm[1]*2,sizeof(double));
    scratch2 = (double *)mxCalloc(dm[0]*dm[1]*2,sizeof(double));
    scratch3 = (double *)mxCalloc(dm[0]*dm[1]*2,sizeof(double));

    cgs2((int *)dm, A, b, rtype, param, tol, nit, x,scratch1,scratch2,scratch3);

    mxFree((void *)scratch3);
    mxFree((void *)scratch2);
    mxFree((void *)scratch1);
}

void fmg2_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int *dm;
    int        cyc=1, nit=1, rtype=2;
    double     *A, *b, *x, *scratch;
    static double param[] = {1.0, 1.0, 1.0, 0.0, 0.0};

    if (nrhs!=3 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfDimensions(prhs[0])!=3) mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetDimensions(prhs[0])[2]!=3)
        mexErrMsgTxt("3rd dimension of 1st arg must be 3.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfDimensions(prhs[1])!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[1]);
    if (dm[2]!=2)
        mexErrMsgTxt("3rd dimension of second arg must be 2.");

    if (mxGetDimensions(prhs[0])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[0])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");

    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Data must be numeric, real, full and double");

    if (mxGetNumberOfElements(prhs[2]) >8)
        mexErrMsgTxt("Third argument should contain rtype, vox1, vox2, param1, param2, param3, ncycles and relax-its.");
    if (mxGetNumberOfElements(prhs[2]) >=1) rtype    = (int)mxGetPr(prhs[2])[0];
    if (mxGetNumberOfElements(prhs[2]) >=2) param[0] = 1/mxGetPr(prhs[2])[1];
    if (mxGetNumberOfElements(prhs[2]) >=3) param[1] = 1/mxGetPr(prhs[2])[2];
    if (mxGetNumberOfElements(prhs[2]) >=4) param[2] = mxGetPr(prhs[2])[3];
    if (mxGetNumberOfElements(prhs[2]) >=5) param[3] = mxGetPr(prhs[2])[4];
    if (mxGetNumberOfElements(prhs[2]) >=6) param[4] = mxGetPr(prhs[2])[5];
    if (mxGetNumberOfElements(prhs[2]) >=7) cyc      = mxGetPr(prhs[2])[6];
    if (mxGetNumberOfElements(prhs[2]) >=8) nit      = (int)mxGetPr(prhs[2])[7];

    plhs[0] = mxCreateNumericArray(3,dm, mxDOUBLE_CLASS, mxREAL);

    A       = mxGetPr(prhs[0]);
    b       = mxGetPr(prhs[1]);
    x       = mxGetPr(plhs[0]);
    scratch = (double *)mxCalloc(fmg2_scratchsize((int *)dm),sizeof(double));
    fmg2((int *)dm, A, b, rtype, param, cyc, nit, x, scratch);
    mxFree((void *)scratch);
}


void vel2mom_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd;
    const int *dm;
    int rtype = 0;
    static double param[] = {1.0, 1.0, 1.0, 0.0, 0.0};

    if (nrhs!=2 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[2]!=2)
        mexErrMsgTxt("3rd dimension must be 2.");

    if (mxGetNumberOfElements(prhs[2]) >6)
        mexErrMsgTxt("Third argument should contain rtype, vox1, vox2, param1, param2, and param3.");
    if (mxGetNumberOfElements(prhs[2]) >=1) rtype    = (int)mxGetPr(prhs[2])[0];
    if (mxGetNumberOfElements(prhs[2]) >=2) param[0] = 1/mxGetPr(prhs[2])[1];
    if (mxGetNumberOfElements(prhs[2]) >=3) param[1] = 1/mxGetPr(prhs[2])[2];
    if (mxGetNumberOfElements(prhs[2]) >=4) param[2] = mxGetPr(prhs[2])[3];
    if (mxGetNumberOfElements(prhs[2]) >=5) param[3] = mxGetPr(prhs[2])[4];
    if (mxGetNumberOfElements(prhs[2]) >=6) param[4] = mxGetPr(prhs[2])[5];

    plhs[0] = mxCreateNumericArray(nd,dm, mxDOUBLE_CLASS, mxREAL);

    if (rtype==1)
        LtLf_me((int *)dm, mxGetPr(prhs[0]), param, mxGetPr(plhs[0]));
    else if (rtype==2)
        LtLf_be((int *)dm, mxGetPr(prhs[0]), param, mxGetPr(plhs[0]));
    else
        LtLf_le((int *)dm, mxGetPr(prhs[0]), param, mxGetPr(plhs[0]));
}


void rsz_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nc[2];
    const int *na;
    double *a, *b, *c;
    if ((nrhs!=2) || (nlhs>1))
        mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
            mexErrMsgTxt("Data must be numeric, real, full and double");

    if (mxGetNumberOfDimensions(prhs[0])!=2) mexErrMsgTxt("Wrong number of dimensions.");
    na = mxGetDimensions(prhs[0]);

    if (mxGetNumberOfElements(prhs[1]) != 2)
    {
        mexErrMsgTxt("Dimensions argument is wrong size.");
    }
    nc[0] = (int)mxGetPr(prhs[1])[0];
    nc[1] = (int)mxGetPr(prhs[1])[1];

    a = mxGetPr(prhs[0]);
    b = (double *)mxCalloc(na[0]*nc[1],sizeof(double));
    plhs[0] = mxCreateNumericArray(2,nc, mxDOUBLE_CLASS, mxREAL);
    c = mxGetPr(plhs[0]);
    resize((int *)na, a, nc, c, b);
    (void)mxFree(b);
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
            fmg2_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"cgs")  || !strcmp(fnc_str,"CGS"))
        {
            mxFree(fnc_str);
            cgs2_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
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
        fmg2_mexFunction(nlhs, plhs, nrhs, prhs);
    }
}

