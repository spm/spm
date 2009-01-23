/* $Id: dartel2.c 2644 2009-01-23 13:01:50Z john $ */
/* (c) John Ashburner */

#include "mex.h"
#include <math.h>
#include "optimizer2d.h"
#include "diffeo2d.h"

static void dartel_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int *dm;
    int        i, k=10, cycles=5, rtype=2, issym=1, nits = 1;
    double     lmreg=0.01, *v, *g, *f, *ov, *scratch, *dj = (double *)0, *ll;
    static int nll[] = {1, 3};
    static double param[] = {1.0, 1.0, 1.0, 0.0, 0.0};
    int dm1[4];

    if ((nrhs!=4 && nrhs!=5) || nlhs>2)
        mexErrMsgTxt("Incorrect usage");

    for(i=0; i<4; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and double");

    if (mxGetNumberOfDimensions(prhs[0])!=3)
        mexErrMsgTxt("Wrong number of dimensions.");
    i = mxGetNumberOfDimensions(prhs[1]);
    if (i!=2 && i!=3)
        mexErrMsgTxt("Wrong number of dimensions.");
    if (mxGetNumberOfDimensions(prhs[2])!=i)
        mexErrMsgTxt("Wrong number of dimensions.");
    dm  = (int *)mxGetDimensions(prhs[0]);
    dm1[0] = mxGetDimensions(prhs[1])[0];
    dm1[1] = mxGetDimensions(prhs[1])[1];
    if (i>2)
        dm1[2] = mxGetDimensions(prhs[1])[2];
    else
        dm1[2] = 1;

    if (mxGetDimensions(prhs[0])[2]!=2)
        mexErrMsgTxt("3rd dimension of 1st arg must be 2.");

    if (dm1[0] != dm[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (dm1[1] != dm[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");

    if (mxGetDimensions(prhs[2])[0] != dm1[0])
        mexErrMsgTxt("Incompatible 1st dimension.");
    if (mxGetDimensions(prhs[2])[1] != dm1[1])
        mexErrMsgTxt("Incompatible 2nd dimension.");

    if (i>2 && mxGetDimensions(prhs[2])[2] != dm1[2])
        mexErrMsgTxt("Incompatible 3rd dimension.");

    if (mxGetNumberOfElements(prhs[3]) >9)
        mexErrMsgTxt("Fourth argument should contain rtype, param1, param2, param3, LMreg, ncycles, nits, K & issym.");
    if (mxGetNumberOfElements(prhs[3]) >=1) rtype    = mxGetPr(prhs[3])[0];
    if (mxGetNumberOfElements(prhs[3]) >=2) param[2] = mxGetPr(prhs[3])[1];
    if (mxGetNumberOfElements(prhs[3]) >=3) param[3] = mxGetPr(prhs[3])[2];
    if (mxGetNumberOfElements(prhs[3]) >=4) param[4] = mxGetPr(prhs[3])[3];
    if (mxGetNumberOfElements(prhs[3]) >=5) lmreg    = mxGetPr(prhs[3])[4];
    if (mxGetNumberOfElements(prhs[3]) >=6) cycles   = mxGetPr(prhs[3])[5];
    if (mxGetNumberOfElements(prhs[3]) >=7) nits     = mxGetPr(prhs[3])[6];
    if (mxGetNumberOfElements(prhs[3]) >=8) k        = mxGetPr(prhs[3])[7];
    if (mxGetNumberOfElements(prhs[3]) >=9) issym    = mxGetPr(prhs[3])[8];

    plhs[0] = mxCreateNumericArray(3, (unsigned int *)dm,  mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, (unsigned int *)nll, mxDOUBLE_CLASS, mxREAL);

    v       = mxGetPr(prhs[0]);
    g       = mxGetPr(prhs[1]);
    f       = mxGetPr(prhs[2]);
    if (nrhs>=5)
        dj  = mxGetPr(prhs[4]);
    ov      = mxGetPr(plhs[0]);
    ll      = mxGetPr(plhs[1]);

    scratch = mxCalloc(dartel_scratchsize((int *)dm1,issym),sizeof(double));
    dartel((int *)dm1, k, v, g, f, dj, rtype, param, lmreg, cycles, nits, issym, ov,ll, scratch);
    mxFree((void *)scratch);
}


static void cgs2_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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
    dm = (int *)mxGetDimensions(prhs[1]);
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

    plhs[0] = mxCreateNumericArray(3, (unsigned int *)dm, mxDOUBLE_CLASS, mxREAL);

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

static void fmg2_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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
    dm = (int *)mxGetDimensions(prhs[1]);
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

    plhs[0] = mxCreateNumericArray(3, (unsigned int *)dm, mxDOUBLE_CLASS, mxREAL);

    A       = mxGetPr(prhs[0]);
    b       = mxGetPr(prhs[1]);
    x       = mxGetPr(plhs[0]);
    scratch = (double *)mxCalloc(fmg2_scratchsize((int *)dm),sizeof(double));
    fmg2((int *)dm, A, b, rtype, param, cyc, nit, x, scratch);
    mxFree((void *)scratch);
}


static void rsz_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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
    na = (int *)mxGetDimensions(prhs[0]);

    if (mxGetNumberOfElements(prhs[1]) != 2)
    {
        mexErrMsgTxt("Dimensions argument is wrong size.");
    }
    nc[0] = (int)mxGetPr(prhs[1])[0];
    nc[1] = (int)mxGetPr(prhs[1])[1];

    a = mxGetPr(prhs[0]);
    b = (double *)mxCalloc(na[0]*nc[1],sizeof(double));
    plhs[0] = mxCreateNumericArray(2, (unsigned int *)nc, mxDOUBLE_CLASS, mxREAL);
    c = mxGetPr(plhs[0]);
    resize((int *)na, a, nc, c, b);
    (void)mxFree(b);
}

static void vel2mom_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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
    dm = (int *)mxGetDimensions(prhs[0]);
    if (dm[2]!=2)
        mexErrMsgTxt("3rd dimension must be 2.");

    if (mxGetNumberOfElements(prhs[1]) >6)
        mexErrMsgTxt("Third argument should contain rtype, vox1, vox2, param1, param2, and param3.");
    if (mxGetNumberOfElements(prhs[1]) >=1) rtype    = (int)mxGetPr(prhs[1])[0];
    if (mxGetNumberOfElements(prhs[1]) >=2) param[0] = 1/mxGetPr(prhs[1])[1];
    if (mxGetNumberOfElements(prhs[1]) >=3) param[1] = 1/mxGetPr(prhs[1])[2];
    if (mxGetNumberOfElements(prhs[1]) >=4) param[2] = mxGetPr(prhs[1])[3];
    if (mxGetNumberOfElements(prhs[1]) >=5) param[3] = mxGetPr(prhs[1])[4];
    if (mxGetNumberOfElements(prhs[1]) >=6) param[4] = mxGetPr(prhs[1])[5];

    plhs[0] = mxCreateNumericArray(nd, (unsigned int *)dm, mxDOUBLE_CLASS, mxREAL);

    if (rtype==1)
        LtLf_me((int *)dm, mxGetPr(prhs[0]), param, mxGetPr(plhs[0]));
    else if (rtype==2)
        LtLf_be((int *)dm, mxGetPr(prhs[0]), param, mxGetPr(plhs[0]));
    else
        LtLf_le((int *)dm, mxGetPr(prhs[0]), param, mxGetPr(plhs[0]));
}

static void comp_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C;
    int nd, m,n, i;
    const int *dm;

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
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and double");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = (int *)mxGetDimensions(prhs[0]);
    m = dm[0];
    n = dm[1];
    if (dm[2]!=2)
        mexErrMsgTxt("3rd dimension must be 2.");

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = (int *)mxGetDimensions(prhs[1]);
    if (dm[0]!=m || dm[1]!=n || dm[2]!=2)
        mexErrMsgTxt("Incompatible dimensions.");

    plhs[0] = mxCreateNumericArray(nd, (unsigned int *)dm, mxDOUBLE_CLASS, mxREAL);

    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(plhs[0]);

    if (nrhs==2)
    {
        (void)composition((int *)dm,B,A,C);
    }
    else if (nrhs==4)
    {
        double *JA, *JB, *JC;
        nd = mxGetNumberOfDimensions(prhs[2]);
        if (nd==2)
        {
            dm = (int *)mxGetDimensions(prhs[2]);
            if (dm[0]!=m || dm[1]!=n)
                mexErrMsgTxt("Incompatible dimensions.");

            nd = mxGetNumberOfDimensions(prhs[3]);
            if (nd!=2) mexErrMsgTxt("Wrong number of dimensions.");
            dm = (int *)mxGetDimensions(prhs[3]);
            if (dm[0]!=m || dm[1]!=n)
                mexErrMsgTxt("Incompatible dimensions.");

            plhs[1] = mxCreateNumericArray(nd, (unsigned int *)dm, mxDOUBLE_CLASS, mxREAL);

            JA = mxGetPr(prhs[2]);
            JB = mxGetPr(prhs[3]);
            JC = mxGetPr(plhs[1]);
            composition_detjac((int *)dm, B, JB, A, JA, C, JC);
        }
        else if (nd==4)
        {
            dm = (int *)mxGetDimensions(prhs[2]);
            if (dm[0]!=m || dm[1]!=n || dm[2]!=2 || dm[3]!=2)
                mexErrMsgTxt("Incompatible dimensions.");

            nd = mxGetNumberOfDimensions(prhs[3]);
            if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
            dm = (int *)mxGetDimensions(prhs[3]);
            if (dm[0]!=m || dm[1]!=n || dm[2]!=2)
                mexErrMsgTxt("Incompatible dimensions.");

            plhs[1] = mxCreateNumericArray(nd,(unsigned int *)dm, mxDOUBLE_CLASS, mxREAL);

            JA = mxGetPr(prhs[2]);
            JB = mxGetPr(prhs[3]);
            JC = mxGetPr(plhs[1]);
            composition_jacobian((int *)dm, B, JB, A, JA, C, JC);
        }
        else
            mexErrMsgTxt("Wrong number of dimensions.");
    }
    unwrap((int *)dm, C);
}

static void exp_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int k=10, nd;
    const int *dm;
    double *v, *t, *t1;

    if (((nrhs != 1) && (nrhs != 2)) || (nlhs>2)) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = (int *)mxGetDimensions(prhs[0]);
    if (dm[2]!=2)
        mexErrMsgTxt("3rd dimension must be 2.");

    if (nrhs>1)
    {
        if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[1]) != 1)
            mexErrMsgTxt("Data must be scalar");
        k   = (int)(mxGetPr(prhs[1])[0]);
    }

    v       = mxGetPr(prhs[0]);

    plhs[0] = mxCreateNumericArray(nd,(unsigned int *)dm, mxDOUBLE_CLASS, mxREAL);
    t       = mxGetPr(plhs[0]);
    t1      = mxCalloc(dm[0]*dm[1]*2,sizeof(double));

    if (nlhs < 2)
    {
        expdef((int *)dm, k, v, t, t1, (double *)0, (double *)0);
    }
    else
    {
        double *J, *J1;
        int dmj[4];
        dmj[0]  = dm[0];
        dmj[1]  = dm[1];
        dmj[2]  = 2;
        dmj[3]  = 2;
        plhs[1] = mxCreateNumericArray(4,(unsigned int*)dmj, mxDOUBLE_CLASS, mxREAL);
        J       = mxGetPr(plhs[1]);
        J1      = mxCalloc(dm[0]*dm[1]*2*2,sizeof(double));
        expdef((int *)dm, k, v, t, t1, J, J1);
        mxFree((void *)J1);
    }
    unwrap((int *)dm, t);
    mxFree((void *)t1);
}

static void expdet_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int k=10, nd;
    const int *dm;
    double *v, *t, *t1;

    if (((nrhs != 1) && (nrhs != 2)) || (nlhs>2)) mexErrMsgTxt("Incorrect usage.");
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = (int *)mxGetDimensions(prhs[0]);
    if (dm[2]!=2)
        mexErrMsgTxt("3rd dimension must be 2.");

    if (nrhs>1)
    {
        if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[1]) != 1)
            mexErrMsgTxt("Data must be scalar");
        k   = (int)(mxGetPr(prhs[1])[0]);
    }

    v       = mxGetPr(prhs[0]);

    plhs[0] = mxCreateNumericArray(nd,(unsigned int *)dm, mxDOUBLE_CLASS, mxREAL);
    t       = mxGetPr(plhs[0]);
    t1      = mxCalloc(dm[0]*dm[1]*2,sizeof(double));

    if (nlhs < 2)
    {
        expdef((int *)dm, k, v, t, t1, (double *)0, (double *)0);
    }
    else
    {
        double *J, *J1;
        int dmj[2];
        dmj[0]  = dm[0];
        dmj[1]  = dm[1];
        plhs[1] = mxCreateNumericArray(2,(unsigned int *)dmj, mxDOUBLE_CLASS, mxREAL);
        J       = mxGetPr(plhs[1]);
        J1      = mxCalloc(dm[0]*dm[1],sizeof(double));
        expdefdet((int *)dm, k, v, t, t1, J, J1);
        mxFree((void *)J1);
    }
    unwrap((int *)dm, t);
    mxFree((void *)t1);
}

static void brc_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C;
    int nd, m,n, i;
    const int *dm;

    if (nrhs == 0) mexErrMsgTxt("Incorrect usage");
    if (nlhs > 1) mexErrMsgTxt("Only 1 output argument required");

    for(i=0; i<nrhs; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and double");

    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = (int *)mxGetDimensions(prhs[0]);
    m = dm[0];
    n = dm[1];
    if (dm[2]!=2)
        mexErrMsgTxt("3rd dimension must be 2.");

    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=3) mexErrMsgTxt("Wrong number of dimensions.");
    dm = (int *)mxGetDimensions(prhs[1]);
    if (dm[0]!=m || dm[1]!=n || dm[2]!=2)
        mexErrMsgTxt("Incompatible dimensions.");

    plhs[0] = mxCreateNumericArray(nd,(unsigned int*)dm, mxDOUBLE_CLASS, mxREAL);

    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(plhs[0]);

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
        if (!strcmp(fnc_str,"dartel") || !strcmp(fnc_str,"DARTEL"))
        {
            mxFree(fnc_str);
            dartel_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"comp"))
        {
            mxFree(fnc_str);
            comp_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"vel2mom") || !strcmp(fnc_str,"mom"))
        {
            mxFree(fnc_str);
            vel2mom_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"Exp") || !strcmp(fnc_str,"exp") || !strcmp(fnc_str,"EXP"))
        {
            mxFree(fnc_str);
            exp_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"Expdet") || !strcmp(fnc_str,"expdet") || !strcmp(fnc_str,"EXPDET"))
        {
            mxFree(fnc_str);
            expdet_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"samp"))
        {
            mxFree(fnc_str);
            /* samp_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]); */
            mexErrMsgTxt("samp is not ready yet.");
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

