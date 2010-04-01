/* $Id: dartel3.c 3805 2010-04-01 15:34:54Z john $ */
/* (c) John Ashburner (2007) */

#include "mex.h"
#include <math.h>
#include "optimizer3d.h"
#include "diffeo3d.h"

static void dartel_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int        i, k=10, cycles=4, its=2, rtype=2, code=0;
     int dm[5];
    double     lmreg0=0.0, *ll;
    float      *v, *g, *f, *jd = (float *)0, *ov, *scratch;
    static double param[6] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    static  int nll[] = {1, 3, 1};
    
    if ((nrhs!=4 && nrhs!=5) || nlhs>2)
        mexErrMsgTxt("Incorrect usage");

    for(i=0; i<3; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
            mexErrMsgTxt("Data must be numeric, real, full and double");

    if (mxGetNumberOfDimensions(prhs[0])!=4) mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetNumberOfDimensions(prhs[1])>4) mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetNumberOfDimensions(prhs[2])!=mxGetNumberOfDimensions(prhs[1])) mexErrMsgTxt("Incompatible number of dimensions.");
    dm[0] = mxGetDimensions(prhs[0])[0];
    dm[1] = mxGetDimensions(prhs[0])[1];
    dm[2] = mxGetDimensions(prhs[0])[2];
    dm[3] = mxGetDimensions(prhs[0])[3];

    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension of 1st arg must be 3.");

    if (mxGetDimensions(prhs[1])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[1])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");
    if (mxGetNumberOfDimensions(prhs[1])>=3 && mxGetDimensions(prhs[1])[2] != dm[2])
        mexErrMsgTxt("Incompatible 3rd dimension.");

    if (mxGetDimensions(prhs[2])[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[2])[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");
    if (mxGetNumberOfDimensions(prhs[2])>=3 && mxGetDimensions(prhs[2])[2] != dm[2])
        mexErrMsgTxt("Incompatible 3rd dimension.");

    if (nrhs>=5)
    {
        if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) || mxIsSparse(prhs[4]) || !mxIsSingle(prhs[4]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
        if (mxGetNumberOfDimensions(prhs[4])!=3) mexErrMsgTxt("Wrong number of dimensions.");
        if (mxGetDimensions(prhs[4])[0] != dm[0])
            mexErrMsgTxt("Incompatible 1st dimension.");
        if (mxGetDimensions(prhs[4])[1] != dm[1])
            mexErrMsgTxt("Incompatible 2nd dimension.");
        if (mxGetDimensions(prhs[4])[2] != dm[2])
            mexErrMsgTxt("Incompatible 3rd dimension.");
        jd = (float *)mxGetPr(prhs[4]);
    }
    if (mxGetNumberOfElements(prhs[3]) >9)
        mexErrMsgTxt("Fourth argument should contain rtype, param1, param2, param3, LMreg, ncycles, nits, nsamps and code.");
    if (mxGetNumberOfElements(prhs[3]) >=1) rtype  = mxGetPr(prhs[3])[0];
    if (mxGetNumberOfElements(prhs[3]) >=2) param[3] = mxGetPr(prhs[3])[1];
    if (mxGetNumberOfElements(prhs[3]) >=3) param[4] = mxGetPr(prhs[3])[2];
    if (mxGetNumberOfElements(prhs[3]) >=4) param[5] = mxGetPr(prhs[3])[3];
    if (mxGetNumberOfElements(prhs[3]) >=5) lmreg0 = mxGetPr(prhs[3])[4];
    if (mxGetNumberOfElements(prhs[3]) >=6) cycles = mxGetPr(prhs[3])[5];
    if (mxGetNumberOfElements(prhs[3]) >=7) its    = mxGetPr(prhs[3])[6];
    if (mxGetNumberOfElements(prhs[3]) >=8) k      = mxGetPr(prhs[3])[7];
    if (mxGetNumberOfElements(prhs[3]) >=9) code   = mxGetPr(prhs[3])[8];

    plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2,nll, mxDOUBLE_CLASS, mxREAL);
    
    v       = (float *)mxGetPr(prhs[0]);
    g       = (float *)mxGetPr(prhs[1]);
    f       = (float *)mxGetPr(prhs[2]);
    ov      = (float *)mxGetPr(plhs[0]);
    ll      = (double*)mxGetPr(plhs[1]);
    
    scratch = (float *)mxCalloc(iteration_scratchsize((int *)dm, code,k),sizeof(float));

    dm[3] = 1;
    if (mxGetNumberOfDimensions(prhs[1])>=4)
        dm[3] = mxGetDimensions(prhs[1])[3];

    iteration(dm, k, v, g, f, jd, rtype, param, lmreg0, cycles, its, code,
              ov, ll, scratch);
    mxFree((void *)scratch);
}


static void cgs3_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const  int *dm;
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

static void fmg3_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const  int *dm;
    int        cyc=1, nit=1, rtype=0;
    float      *A, *b, *x, *scratch;
    static double param[6] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};

    if ((nrhs!=3 && nrhs!=4) || nlhs>1)
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

    if (nrhs>=4)
    {
        int i;
        float *x_orig;
        if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsSingle(prhs[3]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
        if (mxGetNumberOfDimensions(prhs[3])!=4) mexErrMsgTxt("Wrong number of dimensions.");
        if (mxGetDimensions(prhs[3])[3]!=3)
            mexErrMsgTxt("4th dimension of fourth arg must be 3.");
        if (mxGetDimensions(prhs[3])[0] != dm[0])
            mexErrMsgTxt("Incompatible 1st dimension.");
        if (mxGetDimensions(prhs[3])[1] != dm[1])
            mexErrMsgTxt("Incompatible 2nd dimension.");
        if (mxGetDimensions(prhs[3])[2] != dm[2])
            mexErrMsgTxt("Incompatible 3rd dimension.");

        plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);
        x_orig  = (float *)mxGetPr(prhs[3]);
        x       = (float *)mxGetPr(plhs[0]);
        for(i=0; i<dm[0]*dm[1]*dm[2]*3; i++)
            x[i] = x_orig[i];
    }
    else
    {
        plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);
        x       = (float *)mxGetPr(plhs[0]);
    }

    A       = (float *)mxGetPr(prhs[0]);
    b       = (float *)mxGetPr(prhs[1]);
    scratch = (float *)mxCalloc(fmg3_scratchsize((int *)dm),sizeof(float));
    fmg3((int *)dm, A, b, rtype, param, cyc, nit, x, scratch);
    mxFree((void *)scratch);
}

static void fmg3_noa_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const  int *dm;
    int        cyc=1, nit=1, rtype=0;
    float      *b, *x, *scratch;
    static double param[6] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};

    if ((nrhs!=2 && nrhs!=3) || nlhs>1)
        mexErrMsgTxt("Incorrect usage");

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
        mexErrMsgTxt("Data must be numeric, real, full and single");
    if (mxGetNumberOfDimensions(prhs[0])!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension of 1st arg must be 3.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Data must be numeric, real, full and double");

    if (mxGetNumberOfElements(prhs[1]) != 9)
        mexErrMsgTxt("Third argument should contain rtype, vox1, vox2, vox3, param1, param2, param3, ncycles and relax-its.");
    rtype    = (int)(mxGetPr(prhs[1])[0]);
    if (rtype!=0) mexErrMsgTxt("Only does linear elastic energy.");
    param[0] = 1/mxGetPr(prhs[1])[1];
    param[1] = 1/mxGetPr(prhs[1])[2];
    param[2] = 1/mxGetPr(prhs[1])[3];
    param[3] = mxGetPr(prhs[1])[4];
    param[4] = mxGetPr(prhs[1])[5];
    param[5] = mxGetPr(prhs[1])[6];
    cyc      = mxGetPr(prhs[1])[7];
    nit      = (int)(mxGetPr(prhs[1])[8]);

    if (nrhs>=3)
    {
        int i;
        float *x_orig;
        if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsSingle(prhs[2]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
        if (mxGetNumberOfDimensions(prhs[2])!=4) mexErrMsgTxt("Wrong number of dimensions.");
        if (mxGetDimensions(prhs[2])[3]!=3)
            mexErrMsgTxt("4th dimension of third arg must be 3.");
        if (mxGetDimensions(prhs[2])[0] != dm[0])
            mexErrMsgTxt("Incompatible 1st dimension.");
        if (mxGetDimensions(prhs[2])[1] != dm[1])
            mexErrMsgTxt("Incompatible 2nd dimension.");
        if (mxGetDimensions(prhs[2])[2] != dm[2])
            mexErrMsgTxt("Incompatible 3rd dimension.");

        plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);
        x_orig  = (float *)mxGetPr(prhs[2]);
        x       = (float *)mxGetPr(plhs[0]);
        for(i=0; i<dm[0]*dm[1]*dm[2]*3; i++)
            x[i] = x_orig[i];
    }
    else
    {
        plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);
        x       = (float *)mxGetPr(plhs[0]);
    }

    b       = (float *)mxGetPr(prhs[0]);
    scratch = (float *)mxCalloc(fmg3_scratchsize_noa((int *)dm),sizeof(float));
    fmg3_noa((int *)dm, b, rtype, param, cyc, nit, x, scratch);
    mxFree((void *)scratch);
}

static void rsz_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     int nc[3];
    const  int *na;
    float *a, *b, *c;
    if ((nrhs!=2) || (nlhs>1))
        mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
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
    b = (float *)mxCalloc(3*nc[0]*nc[1]+na[0]*nc[1],sizeof(float));
    plhs[0] = mxCreateNumericArray(3,nc, mxSINGLE_CLASS, mxREAL);
    c = (float *)mxGetPr(plhs[0]);
    resize((int *)na, a, nc, c, b);
    (void)mxFree(b);
}

static void vel2mom_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd;
    const  int *dm;
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
        mexErrMsgTxt("Parameters should contain rtype, vox1, vox2, vox3, param1, param2 and param3.");
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

static void comp_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float *A, *B, *C;
    int nd, n0, n1, n2, i;
    const  int *dm;

    if (nrhs == 0) mexErrMsgTxt("Incorrect usage");
    if (nrhs == 2)
    {
        if (nlhs > 1) mexErrMsgTxt("Only 1 output argument required");
    }
    else if (nrhs == 4)
    {
        if (nlhs > 2) mexErrMsgTxt("Only 2 output argument required");
    }
    else
        mexErrMsgTxt("Either 2 or 4 input arguments required");

    for(i=0; i<nrhs; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    n0 = dm[0];
    n1 = dm[1];
    n2 = dm[2];
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension must be 3.");

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[1]);
    if (dm[0]!=n0 || dm[1]!=n1 || dm[2]!=n2 || dm[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");

    plhs[0] = mxCreateNumericArray(nd,dm, mxSINGLE_CLASS, mxREAL);

    A = (float *)mxGetPr(prhs[0]);
    B = (float *)mxGetPr(prhs[1]);
    C = (float *)mxGetPr(plhs[0]);

    if (nrhs==2)
    {
        (void)composition((int *)dm,B,A,C);
    }
    else if (nrhs==4)
    {
        float *JA, *JB, *JC;
        nd = mxGetNumberOfDimensions(prhs[2]);
        if (nd==5)
        {
            dm = mxGetDimensions(prhs[2]);
            if (dm[0]!=n0 || dm[1]!=n1 || dm[2]!=n2 || dm[3]!=3 || dm[4]!=3)
                mexErrMsgTxt("Incompatible dimensions.");

            nd = mxGetNumberOfDimensions(prhs[3]);
            if (nd!=5) mexErrMsgTxt("Wrong number of dimensions.");
            dm = mxGetDimensions(prhs[3]);
            if (dm[0]!=n0 || dm[1]!=n1 || dm[2]!=n2 || dm[3]!=3 || dm[4]!=3)
                mexErrMsgTxt("Incompatible dimensions.");

            plhs[1] = mxCreateNumericArray(nd,dm, mxSINGLE_CLASS, mxREAL);

            JA = (float *)mxGetPr(prhs[2]);
            JB = (float *)mxGetPr(prhs[3]);
            JC = (float *)mxGetPr(plhs[1]);
            composition_jacobian((int *)dm, B, JB, A, JA, C, JC);
        }
        else if (nd==3)
        {
            dm = mxGetDimensions(prhs[2]);
            if (dm[0]!=n0 || dm[1]!=n1 || dm[2]!=n2)
                mexErrMsgTxt("Incompatible dimensions.");

            nd = mxGetNumberOfDimensions(prhs[3]);
            if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
            dm = mxGetDimensions(prhs[3]);
            if (dm[0]!=n0 || dm[1]!=n1 || dm[2]!=n2)
                mexErrMsgTxt("Incompatible dimensions.");

            plhs[1] = mxCreateNumericArray(nd,dm, mxSINGLE_CLASS, mxREAL);

            JA = (float *)mxGetPr(prhs[2]);
            JB = (float *)mxGetPr(prhs[3]);
            JC = (float *)mxGetPr(plhs[1]);
            composition_jacdet((int *)dm, B, JB, A, JA, C, JC);
        }
        else mexErrMsgTxt("Wrong number of dimensions.");
    }
    unwrap((int *)dm, C);
}

static void samp_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float *f, *Y, *wf;
    double buf[1024];
    int nd, i, mm;
     int dmf[4], dmy[4];
    const  int *dmyp;

    if (nrhs == 0) mexErrMsgTxt("Incorrect usage");
    if (nrhs == 2)
    {
        if (nlhs > 1) mexErrMsgTxt("Only 1 output argument required");
    }
    else
        mexErrMsgTxt("Two input arguments required");

    for(i=0; i<nrhs; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    dmf[0] = dmf[1] = dmf[2] = dmf[3] = 1;
    for(i=0; i<nd; i++)
        dmf[i] = mxGetDimensions(prhs[0])[i];

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmyp = mxGetDimensions(prhs[1]);
    if (dmyp[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");

    dmy[0] = dmyp[0];
    dmy[1] = dmyp[1];
    dmy[2] = dmyp[2];
    dmy[3] = dmf[3];
    plhs[0] = mxCreateNumericArray(4,dmy, mxSINGLE_CLASS, mxREAL);

    f = (float *)mxGetPr(prhs[0]);
    Y = (float *)mxGetPr(prhs[1]);
    wf= (float *)mxGetPr(plhs[0]);

    mm = dmy[0]*dmy[1]*dmy[2];
    for (i=0; i<mm; i++)
    {
        int j;
        sampn(dmf, f, dmf[3], mm,
            (double)Y[i]-1.0, (double)Y[mm+i]-1.0, (double)Y[2*mm+i]-1.0,
            buf);
        for(j=0; j<dmf[3]; j++)
            wf[i+mm*j] = buf[j];
    }
}

static void push_mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    float *f, *Y, *so, *po;
    int nd, i, m, n;
     int dmf[4];
     int dmo[4];
    const  int *dmy;

    if ((nrhs != 2) && (nrhs != 3))
        mexErrMsgTxt("Two or three input arguments required");
    if (nlhs  > 2) mexErrMsgTxt("Up to two output arguments required");
    
    for(i=0; i<2; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) ||
             mxIsSparse( prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    dmf[0] = dmf[1] = dmf[2] = dmf[3] = 1;
    for(i=0; i<nd; i++)
        dmf[i] = mxGetDimensions(prhs[0])[i];

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmy = mxGetDimensions(prhs[1]);
    if (dmy[0]!=dmf[0] || dmy[1]!=dmf[1] || dmy[2]!=dmf[2] || dmy[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");
    
    if (nrhs>=3)
    {
        if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) ||
        mxIsSparse( prhs[2]) || !mxIsDouble(prhs[2]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[2])!= 3)
            mexErrMsgTxt("Output dimensions must have three elements");
        dmo[0] = (int)floor(mxGetPr(prhs[2])[0]);
        dmo[1] = (int)floor(mxGetPr(prhs[2])[1]);
        dmo[2] = (int)floor(mxGetPr(prhs[2])[2]);
    }
    else
    {
        dmo[0] = dmf[0];
        dmo[1] = dmf[1];
        dmo[2] = dmf[2];
    }
    dmo[3] = dmf[3];

    plhs[0] = mxCreateNumericArray(4,dmo, mxSINGLE_CLASS, mxREAL);
    f  = (float *)mxGetPr(prhs[0]);
    Y  = (float *)mxGetPr(prhs[1]);
    po = (float *)mxGetPr(plhs[0]);
    if (nlhs>=2)
    {
        plhs[1] = mxCreateNumericArray(3,dmo, mxSINGLE_CLASS, mxREAL);
        so      = (float *)mxGetPr(plhs[1]);
    }
    else
        so      = (float *)0;
    
    m = dmf[0]*dmf[1]*dmf[2];
    n = dmf[3];
    
    push(dmo, m, n, Y, f, po, so);
}

static void pushc_mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    float *f, *Y, *so, *po;
    int nd, i, m, n;
     int dmf[4];
     int dmo[4];
    const  int *dmy;

    if ((nrhs != 2) && (nrhs != 3))
        mexErrMsgTxt("Two or three input arguments required");
    if (nlhs  > 2) mexErrMsgTxt("Up to two output arguments required");

    for(i=0; i<2; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) ||
             mxIsSparse( prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    dmf[0] = dmf[1] = dmf[2] = dmf[3] = 1;
    for(i=0; i<nd; i++)
        dmf[i] = mxGetDimensions(prhs[0])[i];

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmy = mxGetDimensions(prhs[1]);
    if (dmy[0]!=dmf[0] || dmy[1]!=dmf[1] || dmy[2]!=dmf[2] || dmy[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");

    if (nrhs>=3)
    {
        if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) ||
        mxIsSparse( prhs[2]) || !mxIsDouble(prhs[2]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[2])!= 3)
            mexErrMsgTxt("Output dimensions must have three elements");
        dmo[0] = (int)floor(mxGetPr(prhs[2])[0]);
        dmo[1] = (int)floor(mxGetPr(prhs[2])[1]);
        dmo[2] = (int)floor(mxGetPr(prhs[2])[2]);
    }
    else
    {
        dmo[0] = dmf[0];
        dmo[1] = dmf[1];
        dmo[2] = dmf[2];
    }
    dmo[3] = dmf[3];

    plhs[0] = mxCreateNumericArray(4,dmo, mxSINGLE_CLASS, mxREAL);
    f  = (float *)mxGetPr(prhs[0]);
    Y  = (float *)mxGetPr(prhs[1]);
    po = (float *)mxGetPr(plhs[0]);
    if (nlhs>=2)
    {
        plhs[1] = mxCreateNumericArray(3,dmo, mxSINGLE_CLASS, mxREAL);
        so      = (float *)mxGetPr(plhs[1]);
    }
    else
        so      = (float *)0;

    m = dmf[0]*dmf[1]*dmf[2];
    n = dmf[3];

    pushc(dmo, m, n, Y, f, po, so);
}

static void pushc_grads_mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    float *f, *Y, *J, *po;
    int nd, i, m;
    int  dmf[4];
    int  dmo[4];
    const  int *dmy;

    if ((nrhs != 3) && (nrhs != 4))
        mexErrMsgTxt("Two or three input arguments required");
    if (nlhs  > 1) mexErrMsgTxt("Up to one output argument required");

    for(i=0; i<2; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) ||
             mxIsSparse( prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    dmf[0] = dmf[1] = dmf[2] = dmf[3] = 1;
    for(i=0; i<nd; i++)
        dmf[i] = mxGetDimensions(prhs[0])[i];
    if (dmf[3]!=3)
        mexErrMsgTxt("Wrong sized vector field.");

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmy = mxGetDimensions(prhs[1]);
    if (dmy[0]!=dmf[0] || dmy[1]!=dmf[1] || dmy[2]!=dmf[2] || dmy[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");

    nd = mxGetNumberOfDimensions(prhs[2]);
    if (nd!=5) mexErrMsgTxt("Wrong number of dimensions.");
    dmy = mxGetDimensions(prhs[2]);
    if (dmy[0]!=dmf[0] || dmy[1]!=dmf[1] || dmy[2]!=dmf[2] || dmy[3]!=3 || dmy[4]!=3)
        mexErrMsgTxt("Incompatible dimensions.");

    if (nrhs>=4)
    {
        if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) ||
        mxIsSparse( prhs[3]) || !mxIsDouble(prhs[3]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[3])!= 3)
            mexErrMsgTxt("Output dimensions must have three elements");
        dmo[0] = (int)floor(mxGetPr(prhs[3])[0]);
        dmo[1] = (int)floor(mxGetPr(prhs[3])[1]);
        dmo[2] = (int)floor(mxGetPr(prhs[3])[2]);
    }
    else
    {
        dmo[0] = dmf[0];
        dmo[1] = dmf[1];
        dmo[2] = dmf[2];
    }
    dmo[3] = dmf[3];

    plhs[0] = mxCreateNumericArray(4,dmo, mxSINGLE_CLASS, mxREAL);
    f  = (float *)mxGetPr(prhs[0]);
    Y  = (float *)mxGetPr(prhs[1]);
    J  = (float *)mxGetPr(prhs[2]);
    po = (float *)mxGetPr(plhs[0]);

    m = dmf[0]*dmf[1]*dmf[2];
    pushc_grads(dmo, m, Y, J, f, po);
}


static void exp_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int k=6, nd;
    const  int *dm;
    float *v, *t, *t1;
    double sc = 1.0;
    int flg = 0;

    if (((nrhs != 1) && (nrhs != 2)) || (nlhs>2)) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension must be 3.");

    if (nrhs>1)
    {
        if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[1]) > 3)
            mexErrMsgTxt("Params must contain one to three elements");
        if (mxGetNumberOfElements(prhs[1]) >= 1) k   =   (int)(mxGetPr(prhs[1])[0]);
        if (mxGetNumberOfElements(prhs[1]) >= 2) sc  = (float)(mxGetPr(prhs[1])[1]);
        if (mxGetNumberOfElements(prhs[1]) >= 3) flg =   (int)(mxGetPr(prhs[1])[2]);
    }

    v       = (float *)mxGetPr(prhs[0]);

    plhs[0] = mxCreateNumericArray(nd,dm, mxSINGLE_CLASS, mxREAL);
    t       = (float *)mxGetPr(plhs[0]);
    t1      = mxCalloc(dm[0]*dm[1]*dm[2]*3,sizeof(float));

    if (nlhs < 2)
    {
        expdef((int *)dm, k, sc, v, t, t1, (float *)0, (float *)0);
    }
    else
    {
        float *J, *J1;
         int dmj[5];
        dmj[0]  = dm[0];
        dmj[1]  = dm[1];
        dmj[2]  = dm[2];
        if (flg==0)
        {
            dmj[3]  = 3;
            dmj[4]  = 3;
            plhs[1] = mxCreateNumericArray(5,dmj, mxSINGLE_CLASS, mxREAL);
            J       = (float *)mxGetPr(plhs[1]);
            J1      = mxCalloc(dm[0]*dm[1]*dm[2]*3*3,sizeof(float));
            expdef((int *)dm, k, sc, v, t, t1, J, J1);
        }
        else
        {
            plhs[1] = mxCreateNumericArray(3,dmj, mxSINGLE_CLASS, mxREAL);
            J       = (float *)mxGetPr(plhs[1]);
            J1      = mxCalloc(dm[0]*dm[1]*dm[2],sizeof(float));
            expdefdet((int *)dm, k, sc, v, t, t1, J, J1);
        }
        mxFree((void *)J1);
    }
    unwrap((int *)dm, t);
    mxFree((void *)t1);
}

static void smalldef_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd;
    const  int *dm;
    float *v, *t;
    double sc = 1.0;

    if (((nrhs != 1) && (nrhs != 2)) || (nlhs>2)) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension must be 3.");

    if (nrhs>1)
    {
        if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[1]) > 1)
            mexErrMsgTxt("Params must contain one element");
        if (mxGetNumberOfElements(prhs[1]) >= 1) sc  = (float)(mxGetPr(prhs[1])[0]);
    }

    v       = (float *)mxGetPr(prhs[0]);

    plhs[0] = mxCreateNumericArray(nd,dm, mxSINGLE_CLASS, mxREAL);
    t       = (float *)mxGetPr(plhs[0]);

    if (nlhs < 2)
    {
        smalldef((int *)dm, sc, v, t);
    }
    else
    {
        float *J, *J1;
        int dmj[5];
        dmj[0]  = dm[0];
        dmj[1]  = dm[1];
        dmj[2]  = dm[2];
        dmj[3]  = 3;
        dmj[4]  = 3;
        plhs[1] = mxCreateNumericArray(5,dmj, mxSINGLE_CLASS, mxREAL);
        J       = (float *)mxGetPr(plhs[1]);
        smalldef_jac1((int *)dm, sc, v, t, J);
    }
}

static void det_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd;
    const int *dm;

    if ((nrhs != 1) || (nlhs>1)) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=5) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[3]!=3) mexErrMsgTxt("4th dimension must be 3.");
    if (dm[4]!=3) mexErrMsgTxt("5th dimension must be 3.");

    plhs[0] = mxCreateNumericArray(3,dm, mxSINGLE_CLASS, mxREAL);
    determinant((int *)dm,(float *)mxGetPr(prhs[0]),(float *)mxGetPr(plhs[0]));
}

static void minmax_div_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nd;
    const int *dm;
    static  int nout[] = {1, 2, 1};

    if ((nrhs != 1) || (nlhs>1)) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[3]!=3) mexErrMsgTxt("4th dimension must be 3.");

    plhs[0] = mxCreateNumericArray(2,nout, mxDOUBLE_CLASS, mxREAL);
    minmax_div((int *)dm,(float *)mxGetPr(prhs[0]),(double *)mxGetPr(plhs[0]));
}

static void brc_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float *A, *B, *C;
    int nd, i;
    const  int *dm, *dm1;

    if (nrhs == 0) mexErrMsgTxt("Incorrect usage");
    if (nrhs != 2) mexErrMsgTxt("Incorrect number of input arguments");
    if (nlhs > 1) mexErrMsgTxt("Only 1 output argument required");

    for(i=0; i<nrhs; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and double");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm = mxGetDimensions(prhs[0]);
    if (dm[3]!=3)
        mexErrMsgTxt("4th dimension must be 3.");

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dm1 = mxGetDimensions(prhs[1]);
    if (dm[0]!=dm1[0] || dm[1]!=dm1[1] || dm[2]!=dm1[2] || dm[3]!=dm1[3])
        mexErrMsgTxt("Incompatible dimensions.");

    plhs[0] = mxCreateNumericArray(nd,dm, mxSINGLE_CLASS, mxREAL);

    A = (float *)mxGetPr(prhs[0]);
    B = (float *)mxGetPr(prhs[1]);
    C = (float *)mxGetPr(plhs[0]);

    (void)bracket((int *)dm,A,B,C);
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
        if (!strcmp(fnc_str,"dartel") || !strcmp(fnc_str,"DARTEL") )
        {
            mxFree(fnc_str);
            dartel_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"comp"))
        {
            mxFree(fnc_str);
            comp_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"vel2mom"))
        {
            mxFree(fnc_str);
            vel2mom_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"Exp"))
        {
            mxFree(fnc_str);
            exp_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"smalldef"))
        {
            mxFree(fnc_str);
            smalldef_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"samp"))
        {
            mxFree(fnc_str);
            samp_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"push"))
        {
            mxFree(fnc_str);
            push_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"pushc"))
        {
            mxFree(fnc_str);
            pushc_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"pushg"))
        {
            mxFree(fnc_str);
            pushc_grads_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"det"))
        {
            mxFree(fnc_str);
            det_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"divrange"))
        {
            mxFree(fnc_str);
            minmax_div_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"fmg")  || !strcmp(fnc_str,"FMG"))
        {
            mxFree(fnc_str);
            fmg3_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"mom2vel"))
        {
            mxFree(fnc_str);
            fmg3_noa_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"cgs")  || !strcmp(fnc_str,"CGS"))
        {
            mxFree(fnc_str);
            cgs3_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"rsz")  || !strcmp(fnc_str,"resize"))
        {
            mxFree(fnc_str);
            rsz_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"brc")  || !strcmp(fnc_str,"bracket"))
        {
            mxFree(fnc_str);
            brc_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else
        {
            mxFree(fnc_str);
            mexErrMsgTxt("Option not recognised.");
        }
    }
    else
    {
        dartel_mexFunction(nlhs, plhs, nrhs, prhs);
    }
}

