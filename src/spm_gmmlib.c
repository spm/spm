/*
 * John Ashburner, Mikael Brudfors & Yael Balbastre
 * Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging
 */

#include <math.h>
#include <string.h>
#include "mex.h"
#include "gmmlib.h"
/* #include "spm_openmp.h" for later??*/

/*
 * r = resp(m,b,W,nu,gam,lkp,mu,f,E,skip)
 */
static mwSize copy_dims(const mxArray *prhs, mwSize n[])
{
    mwSize i, nd;
    const mwSize *dm;
    nd = mxGetNumberOfDimensions(prhs);
    if (nd>5) mexErrMsgTxt("Too many dimensions.");
    dm = mxGetDimensions(prhs);
    for(i=0; i<nd; i++) n[i] = dm[i];
    for(i=nd; i<5; i++) n[i] = 1;
    return nd;
}

static void parse_rhs(int nrhs, const mxArray *prhs[], mwSize *Kp, double **mp, double **bp, double **Wp, double **nup, double **gamp,
                      mwSize **lkpp, mwSize *nm, float **mup, mwSize *nf, float **mfp, float **vfp, mwSize *skip,
                      unsigned char **labelp, double **lnPp)
{
    mwSize nl[5], nd, i, k, P, K1;
    mwSignedIndex *lkp0;

    if (nrhs<9) mexErrMsgTxt("Incorrect usage");

    for(i=0; i<=4; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
            mexErrMsgTxt("GMM parameters must be numeric, real, full and double");

    /* m */
    nd  = copy_dims(prhs[0],nl);
    if (nd>2) mexErrMsgTxt("Wrong number of dimensions (m).");
    P   = nl[0];
    *Kp = nl[1];
    *mp = mxGetPr(prhs[0]);

    /* b */
    nd  = copy_dims(prhs[1],nl);
    if (nd>2) mexErrMsgTxt("Wrong number of dimensions (b).");
    if (nl[0]!=1 || nl[1]!=*Kp) mexErrMsgTxt("Incompatible dimensions (b).");
    *bp  = mxGetPr(prhs[1]);

    /* W */
    nd  = copy_dims(prhs[2],nl);
    if (nd>3) mexErrMsgTxt("Wrong number of dimensions (W).");
    if (nl[0]!=P || nl[1]!=P || nl[2]!=*Kp) mexErrMsgTxt("Incompatible dimensions (W).");
    *Wp  = mxGetPr(prhs[2]);

    /* nu */
    nd  = copy_dims(prhs[3],nl);
    if (nd>2) mexErrMsgTxt("Wrong number of dimensions (nu).");
    if (nl[0]!=1 || nl[1]!=*Kp) mexErrMsgTxt("Incompatible dimensions (nu).");
    *nup = mxGetPr(prhs[3]);

    /* gam */
    nd  = copy_dims(prhs[4],nl);
    if (nd>2) mexErrMsgTxt("Wrong number of dimensions (gam).");
    if (nl[0]!=1 || nl[1]!=*Kp) mexErrMsgTxt("Incompatible dimensions (gam).");
    *gamp = mxGetPr(prhs[4]);

    for(i=0; i<*Kp; i++)
    {
        if ( (*nup)[i]+1.0 <= (double)P) mexErrMsgTxt("Bad nu value.");
        if (  (*bp)[i]     <= 0.0)       mexErrMsgTxt("Bad b value.");
        if ((*gamp)[i]     <= 0.0)       mexErrMsgTxt("Bad gam value.");
    }

    /* lkp */
    if (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5]) || mxIsSparse(prhs[5]) || !mxIsUint64(prhs[5]))
        mexErrMsgTxt("Lookup data must be numeric, real, full and UInt64.");
    nd  = copy_dims(prhs[5],nl);
    if (nd>2) mexErrMsgTxt("Wrong number of dimensions (lkp).");
    if (nl[0]!=1 || nl[1]!=*Kp) mexErrMsgTxt("Incompatible dimensions (lkp).");
    lkp0 = (mwSignedIndex *)mxGetPr(prhs[5]);


    for(i=6; i<=8; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Image data must be numeric, real, full and single.");

    /* mu */
    nd   = copy_dims(prhs[6],nm);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions (mu).");
    *mup = (float *)mxGetPr(prhs[6]);
    K1   = nm[3];

    for(k=0; k<*Kp; k++)
        if (lkp0[k]<1 || lkp0[k]>K1)
            mexErrMsgTxt("Lookup data out of range.");

    /* mf */
    nd  = copy_dims(prhs[7],nf);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions (mf).");
    if (nf[3]!=P) mexErrMsgTxt("Incompatible dimensions (mf).");
    *mfp = (float *)mxGetPr(prhs[7]);

    /* vf */
    nd  = copy_dims(prhs[8],nl);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions (vf).");
    if (nl[0]!=nf[0] || nl[1]!=nf[1] || nl[2]!=nf[2]) mexErrMsgTxt("Incompatible dimensions (vf).");
    if (nl[3]!=P) mexErrMsgTxt("Incompatible dimensions (vf).");
    *vfp = (float *)mxGetPr(prhs[8]);

    /* skip */
    skip[0] = skip[1] = skip[2] = 1;
    if (nrhs>=10)
    {
        mwSize ds[5], i;
        mwSignedIndex *ptr;
        if (!mxIsNumeric(prhs[9]) || mxIsComplex(prhs[9]) ||
              mxIsSparse(prhs[9]) || !mxIsUint64(prhs[9]))
            mexErrMsgTxt("Skip data must be numeric, real, full and UInt64.");
        nd  = copy_dims(prhs[9],ds);
        if (nd>2) mexErrMsgTxt("Wrong number of dimensions (skip).");
        if (ds[0]*ds[1] > 3) mexErrMsgTxt("Wrong number of elements (skip).");
        ptr = (mwSignedIndex *)mxGetPr(prhs[9]);
        for(i=0; i<ds[0]*ds[1]; i++)
        {
            if(ptr[i] > 0) skip[i] = ptr[i];
        }
    }

    /* labels */
    *labelp = NULL;
    *lnPp   = NULL;
    if (nrhs>=12)
    {
        mwSize ds[5];
        nd  = copy_dims(prhs[10],ds);
        if (ds[0]*ds[1]*ds[2]*ds[3] > 0)
        {
            if (!mxIsNumeric(prhs[10]) || mxIsComplex(prhs[10]) || mxIsSparse(prhs[10]) || !mxIsUint8(prhs[10]))
                mexErrMsgTxt("Label data must be numeric, real, full and UInt8.");
            if (ds[0]!=nm[0] || ds[1]!=nm[1] || ds[2]!=nm[2] || ds[3]!=1 || ds[4]!=1)
                mexErrMsgTxt("Incompatible dimensions (label).");
            *labelp = (unsigned char *)mxGetPr(prhs[10]);

            if (!mxIsNumeric(prhs[11]) || mxIsComplex(prhs[11]) || mxIsSparse(prhs[11]) || !mxIsDouble(prhs[11]))
                mexErrMsgTxt("Label parameters must be numeric, real, full and double");
            nd  = copy_dims(prhs[11],ds);
            if (nd>2 || ds[0]!=256 || ds[1]!=*Kp)
                mexErrMsgTxt("Incompatible dimensions (label parameters).");
            *lnPp = (double *)mxGetPr(prhs[11]);
         }
    }

    /* lkp - allocated memory needs freeing later on */
    *lkpp = (mwSize *)mxCalloc(sizeof(mwSize),*Kp);
    for(k=0; k<*Kp; k++) (*lkpp)[k] = lkp0[k]-1;
}

static void resp_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    static mwSize scalar_dim[2] = {1,1};
    mwSize K, k, nf[5], nm[5], nr[5], skip[3];
    double *m, *b, *W, *nu, *gam, *lnP = NULL, *ll;
    float *mu, *mf, *vf, *r;
    mwSize *lkp;
    unsigned char *label = NULL;

    if (nrhs<9 || nrhs>12 || nlhs>2) mexErrMsgTxt("Incorrect usage");

    parse_rhs(nrhs, prhs, &K, &m, &b, &W, &nu, &gam, &lkp, nm, &mu, nf, &mf, &vf, skip, &label, &lnP);
    if (nf[0]!=nm[0] || nf[1]!=nm[1] || nf[2]!=nm[2]) mexErrMsgTxt("Incompatible dimensions (mf).");

    /* r */
    nr[0]   = nf[0];
    nr[1]   = nf[1];
    nr[2]   = nf[2];
    nr[3]   = nm[3]-1;
    plhs[0] = mxCreateNumericArray(4,nr, mxSINGLE_CLASS, mxREAL);
    r       = (float *)mxGetPr(plhs[0]);

    /* ll */
    plhs[1] = mxCreateNumericArray(2,scalar_dim, mxDOUBLE_CLASS, mxREAL);
    ll      = (double *)mxGetPr(plhs[1]);

    ll[0] = call_responsibilities(nf,skip,mf,vf, label, K,m,b,W,nu,gam, lnP, nm[3],lkp,mu, r);
    mxFree(lkp);
}

static void moments_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize K, K1, k, nf[5], nm[5], skip[3], *lkp, dm0, dm1, dm2;
    double *m, *b, *W, *nu, *gam, *lnP = NULL, *s0, *s1, *s2, *ll, *H;
    float *mu, *mf, *vf, *r;
    unsigned char *label = NULL;

    if (nrhs<9 || nrhs>12 || nlhs>5) mexErrMsgTxt("Incorrect usage");

    parse_rhs(nrhs, prhs, &K, &m, &b, &W, &nu, &gam, &lkp, nm, &mu, nf, &mf, &vf, skip, &label, &lnP);
    space_needed(nf[3], K, &dm0, &dm1, &dm2);
    plhs[0] = mxCreateDoubleMatrix(dm0,1, mxREAL); s0 = (double *)mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(dm1,1, mxREAL); s1 = (double *)mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(dm2,1, mxREAL); s2 = (double *)mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(  1,1, mxREAL); ll = (double *)mxGetPr(plhs[3]);
    if (label!=NULL)
    {
        plhs[4] = mxCreateDoubleMatrix(256, K, mxREAL);
        H       = (double *)mxGetPr(plhs[4]);
    }
    else
    {
        plhs[4] = mxCreateDoubleMatrix(0, 0, mxREAL);
        H       = NULL;
    }
    ll[0]   = call_suffstats_missing(nf,mf,vf, label, K,m,b,W,nu,gam, lnP, nm,skip,lkp,mu, s0,s1,s2, H);
    mxFree(lkp);
}

static void inugrads_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize nf[5];
    mwSize K, K1, k, nm[5], dc[5], skip[3], *lkp, dm0, dm1, dm2, nd;
    mwSize c;
    double *m, *b, *W, *nu, *gam, *lnP = NULL, *ll;
    float *mu, *mf, *vf, *g1, *g2;
    unsigned char *label = NULL;

    if (nrhs>13 || nlhs>3) mexErrMsgTxt("Incorrect usage");
    parse_rhs(nrhs, prhs, &K, &m, &b, &W, &nu, &gam, &lkp, nm, &mu, nf, &mf, &vf, skip, &label, &lnP);

    /* c */
    if (!mxIsNumeric(prhs[12]) || mxIsComplex(prhs[12]) ||
          mxIsSparse(prhs[12]) || !mxIsUint64(prhs[12]))
        mexErrMsgTxt("Index must be numeric, real, full and UInt64.");
    nd  = copy_dims(prhs[12],dc);
    if (nd>2 || dc[0]!=1 || dc[1]!=1) mexErrMsgTxt("Index not a scalar.");
    c = ((mwSize *)mxGetPr(prhs[12]))[0] - 1;
    if (c<0 || c>=(mwSize)nf[3]) mexErrMsgTxt("Index out of range.");

    plhs[0] = mxCreateNumericArray(3,nf, mxSINGLE_CLASS, mxREAL); g1 = (float *)mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(3,nf, mxSINGLE_CLASS, mxREAL); g2 = (float *)mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    ll      = (double *)mxGetPr(plhs[2]);

    ll[0]   = call_INUgrads(nf,mf,vf, label, K,m,b,W,nu,gam, lnP, nm,skip,lkp,mu, (mwSize)c, g1,g2);

    mxFree(lkp);
}

static void infer_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize nf[5];
    mwSize K, K1, k, nm[5], dc[5], skip[3], *lkp, dm0, dm1, dm2, nd;
    double *m, *b, *W, *nu, *gam, *lnP = NULL, *sts;
    float *mu, *mf, *vf, *mf1;
    unsigned char *label = NULL;

    if (nrhs>13 || nlhs>3) mexErrMsgTxt("Incorrect usage");
    parse_rhs(nrhs, prhs, &K, &m, &b, &W, &nu, &gam, &lkp, nm, &mu, nf, &mf, &vf, skip, &label, &lnP);

    plhs[0] = mxCreateNumericArray(4,nf, mxSINGLE_CLASS, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mf1     = (float  *)mxGetPr(plhs[0]);
    sts     = (double *)mxGetPr(plhs[1]);

    sts[0]  = (double)call_fill_missing(nf,mf,vf, label, K, m,b,W,nu,gam, lnP, nm,skip,lkp,mu, mf1);

    mxFree(lkp);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* spm_set_num_threads(spm_get_num_threads()); for later?? */

    if ((nrhs>=1) && mxIsChar(prhs[0]))
    {
        char *fnc_str = mxArrayToString(prhs[0]);

        if (!strcmp(fnc_str,"resp"))
        {
            mxFree(fnc_str);
            resp_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"moments"))
        {
            mxFree(fnc_str);
            moments_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"inugrads"))
        {
            mxFree(fnc_str);
            inugrads_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"infer"))
        {
            mxFree(fnc_str);
            infer_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else
        {
            mxFree(fnc_str);
            mexErrMsgTxt("Option not recognised.");
        }
    }
    else
    {
        resp_mexFunction(nlhs, plhs, nrhs, prhs);
    }
}
