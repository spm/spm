/* $Id: optimizer3d_mex.c 964 2007-10-19 16:35:34Z john $ */
/* (c) John Ashburner (2007) */

#include "mex.h"
#include <math.h>
#include "optimizer3d.h"

void cgs3_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int *dm;
    int        nit=1, rtype=0;
    double     tol=1e-10;
    float      *A, *b, *x, *scratch1, *scratch2, *scratch3;
    static double param[6] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};

    if (nrhs!=3 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and single");
    if (mxGetNumberOfDimensions(prhs[0])!=4) mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetDimensions(prhs[0])[3]!=6)
        mexErrMsgTxt("4th dimension of 1st arg must be 6.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsSingle(prhs[1]))
        mexErrMsgTxt("Data must be numeric, real, full and single");
    if (mxGetNumberOfDimensions(prhs[1])!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[1]);
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension of second arg must be 3.");

    if (mxGetDimensions(prhs[0])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[0])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");
    if (mxGetDimensions(prhs[0])[2] != dm[1])
        mexErrMsgTxt("Incompatible 3rd dimension.");

    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Data must be numeric, real, full and double");
    if (mxGetNumberOfElements(prhs[2]) != 9)
        mexErrMsgTxt("Third argument should contain rtype, vox1, vox2, vox3, param1, param2, param3, tol and nit.");

    rtype    = (int)(mxGetPr(prhs[2])[0]);
    param[0] = 1/mxGetPr(prhs[2])[1];
    param[1] = 1/mxGetPr(prhs[2])[2];
    param[2] = 1/mxGetPr(prhs[2])[3];
    param[3] = mxGetPr(prhs[2])[4];
    param[4] = mxGetPr(prhs[2])[5];
    param[5] = mxGetPr(prhs[2])[6];
    tol      = mxGetPr(prhs[2])[7];
    nit      = (int)(mxGetPr(prhs[2])[8]);

    plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);

    A       = (float *)mxGetPr(prhs[0]);
    b       = (float *)mxGetPr(prhs[1]);
    x       = (float *)mxGetPr(plhs[0]);

    scratch1 = (float *)mxCalloc(dm[0]*dm[1]*dm[2]*3,sizeof(float));
    scratch2 = (float *)mxCalloc(dm[0]*dm[1]*dm[2]*3,sizeof(float));
    scratch3 = (float *)mxCalloc(dm[0]*dm[1]*dm[2]*3,sizeof(float));

    cgs3((int *)dm, A, b, rtype, param, tol, nit, x,scratch1,scratch2,scratch3);

    mxFree((void *)scratch3);
    mxFree((void *)scratch2);
    mxFree((void *)scratch1);
}

void fmg3_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int *dm;
    int        cyc=1, nit=1, rtype=0;
    float      *A, *b, *x, *scratch;
    static double param[6] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};

    if (nrhs!=3 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and single");
    if (mxGetNumberOfDimensions(prhs[0])!=4) mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetDimensions(prhs[0])[3]!=6)
        mexErrMsgTxt("4th dimension of 1st arg must be 6.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsSingle(prhs[1]))
        mexErrMsgTxt("Data must be numeric, real, full and single");
    if (mxGetNumberOfDimensions(prhs[1])!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[1]);
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension of second arg must be 3.");

    if (mxGetDimensions(prhs[0])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[0])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");
    if (mxGetDimensions(prhs[0])[2] != dm[2])
        mexErrMsgTxt("Incompatible 3rd dimension.");

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

    plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);

    A       = (float *)mxGetPr(prhs[0]);
    b       = (float *)mxGetPr(prhs[1]);
    x       = (float *)mxGetPr(plhs[0]);
    scratch = (float *)mxCalloc(fmg3_scratchsize((int *)dm),sizeof(float));
    fmg3((int *)dm, A, b, rtype, param, cyc, nit, x, scratch);
    mxFree((void *)scratch);
}

void vel2mom_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd;
    const int *dm;
    int rtype = 0;
    static double param[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};

    if (nrhs!=2 || nlhs>1)
        mexErrMsgTxt("Incorrect usage");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and single");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension must be 3.");

    if (mxGetNumberOfElements(prhs[1]) != 7)
        mexErrMsgTxt("Third argument should contain rtype, vox1, vox2, vox3, param1, param2 and param3.");
    rtype    = (int)(mxGetPr(prhs[1])[0]);
    param[0] = 1/mxGetPr(prhs[1])[1];
    param[1] = 1/mxGetPr(prhs[1])[2];
    param[2] = 1/mxGetPr(prhs[1])[3];
    param[3] = mxGetPr(prhs[1])[4];
    param[4] = mxGetPr(prhs[1])[5];
    param[5] = mxGetPr(prhs[1])[6];

    plhs[0] = mxCreateNumericArray(nd,dm, mxSINGLE_CLASS, mxREAL);

    if (rtype==1)
        LtLf_me((int *)dm, (float *)mxGetPr(prhs[0]), param, (float *)mxGetPr(plhs[0]));
    else if (rtype==2)
        LtLf_be((int *)dm, (float *)mxGetPr(prhs[0]), param, (float *)mxGetPr(plhs[0]));
    else
        LtLf_le((int *)dm, (float *)mxGetPr(prhs[0]), param, (float *)mxGetPr(plhs[0]));
}

void rsz_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nc[3];
    const int *na;
    float *a, *b, *c;
    if ((nrhs!=2) || (nlhs>1))
        mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and float");
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

    a = (float *)mxGetPr(prhs[0]);
    b = (float *)mxCalloc(na[0]*nc[1]+3*nc[0]*nc[1],sizeof(float));
    plhs[0] = mxCreateNumericArray(3,nc, mxSINGLE_CLASS, mxREAL);
    c = (float *)mxGetPr(plhs[0]);
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
            fmg3_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"cgs")  || !strcmp(fnc_str,"CGS"))
        {
            mxFree(fnc_str);
            cgs3_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
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
        fmg3_mexFunction(nlhs, plhs, nrhs, prhs);
    }
}

