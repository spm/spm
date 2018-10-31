/* $Id: spm_field.c 7464 2018-10-31 16:57:27Z john $ */
/* (c) John Ashburner (2007) */

#include "mex.h"
#include <math.h>
#include "shoot_optimN.h"
#include "shoot_boundary.h"

static void boundary_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ((nlhs<=1) && (nrhs==0))
    {
        mwSize nout[] = {1,1,1};
        plhs[0] = mxCreateNumericArray(2,nout, mxDOUBLE_CLASS, mxREAL);
        mxGetPr(plhs[0])[0] = get_bound();
    }
    else if ((nrhs==1) && (nlhs==0))
    {
        if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        set_bound(mxGetPr(prhs[0])[0]);
    }
}

static void fmg_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd, i;
    mwSize dm[4];
    int   cyc=2, nit=2, nparam;
    float *A, *b, *x;
    static double param[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    double scal[256], t[3];

    if ((nrhs!=2 && nrhs!=3 && nrhs!=4) || nlhs>1)
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

    if (nrhs>=3)
    { 
        if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
            mexErrMsgTxt("Third argument must be numeric, real, full and double");
        nparam = mxGetNumberOfElements(prhs[2]);
    }
    else
        nparam = 0;

    if (nparam > 8)
        mexErrMsgTxt("Third argument should contain vox1, vox2, vox3, param1, param2, param3, ncycles and relax-its.");
    if (nparam>=1) param[0] = 1/mxGetPr(prhs[2])[0]; else param[0] = 1.0;
    if (nparam>=2) param[1] = 1/mxGetPr(prhs[2])[1]; else param[1] = 1.0;
    if (nparam>=3) param[2] = 1/mxGetPr(prhs[2])[2]; else param[2] = 1.0;
    if (nparam>=4) param[3] = mxGetPr(prhs[2])[3];   else param[3] = 0.0;
    if (nparam>=5) param[4] = mxGetPr(prhs[2])[4];   else param[4] = 0.0;
    if (nparam>=6) param[5] = mxGetPr(prhs[2])[5];   else param[5] = 0.0;
    if (nparam>=7) cyc      = (int)(mxGetPr(prhs[2])[6]); else cyc = 2;
    if (nparam>=8) nit      = (int)(mxGetPr(prhs[2])[7]); else nit = 2;

    /* Penalise absolute displacements slightly in case supplied Hessian is too small.
       Extra penalty based on value in centre of difference operator, scaled by some
       slightly arbitrary amount.

       On second thoughts - I need to think about this one a bit more as it uses too
       much regularisation if a decent Hessian is supplied.
    */    
    t[0]     = param[0]*param[0];
    t[1]     = param[1]*param[1];
    t[2]     = param[2]*param[2];

    param[3]+= (param[5]*(6*(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]) + 8*(t[0]*t[1]+t[0]*t[2]+t[1]*t[2]))+param[4]*2*(t[0]+t[1]+t[2]))*4e-7;


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
    plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);

    A       = (float *)mxGetPr(prhs[0]);
    b       = (float *)mxGetPr(prhs[1]);
    x       = (float *)mxGetPr(plhs[0]);
    if (param[4]>0.0 || param[5]>0.0)
    {
        float *scratch;
        scratch = (float *)mxCalloc(fmg_scratchsize(dm),sizeof(float));
        fmg(dm, A, b, param, scal, cyc, nit, x, scratch);
        mxFree((void *)scratch);
    }
    else
        solve(dm, A, b, param, scal, x);

}

static void vel2mom_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd, i;
    mwSize dm[4];
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

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Parameters must be numeric, real, full and double");
    if (mxGetNumberOfElements(prhs[1]) != 6)
        mexErrMsgTxt("Parameters should contain vox1, vox2, vox3, param1, param2 and param3.");
    param[0] = 1/mxGetPr(prhs[1])[0];
    param[1] = 1/mxGetPr(prhs[1])[1];
    param[2] = 1/mxGetPr(prhs[1])[2];
    param[3] = mxGetPr(prhs[1])[3];
    param[4] = mxGetPr(prhs[1])[4];
    param[5] = mxGetPr(prhs[1])[5];

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
 
    plhs[0] = mxCreateNumericArray(nd, dm, mxSINGLE_CLASS, mxREAL);

    LtLf(dm, (float *)mxGetPr(prhs[0]), param, scal, (float *)mxGetPr(plhs[0]));
}

#include<string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    set_bound(get_bound());
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
        else if (!strcmp(fnc_str,"boundary")  || !strcmp(fnc_str,"bound"))
        {
            mxFree(fnc_str);
            boundary_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
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

