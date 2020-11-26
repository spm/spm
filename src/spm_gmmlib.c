/* 
 * Copyright (c) 2020 Wellcome Centre for Human Neuroimaging
 * John Ashburner, Mikael Brudfors & Yael Balbastre
 * $Id: spm_gmmlib.c 8021 2020-11-26 15:47:56Z john $
 *
 */

#include <math.h>
#include <string.h>
#include "mex.h"
#include "gmmlib.h"
#include "spm_openmp.h"

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

static void parse_rhs(int nrhs, const mxArray *prhs[], size_t *Kp, double **mp, double **bp, double **Wp, double **nup, double **gamp,
                      size_t **lkpp, size_t *nm, float **mup, size_t *nf, float **mfp, float **vfp, size_t *skip)
{
    size_t nl[5], nd, i, k, P, K1;
    unsigned long long *lkp0;

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

    /* lkp */
    if (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5]) || mxIsSparse(prhs[5]) || !mxIsUint64(prhs[5]))
        mexErrMsgTxt("Lookup data must be numeric, real, full and UInt64.");
    nd  = copy_dims(prhs[5],nl);
    if (nd>2) mexErrMsgTxt("Wrong number of dimensions (lkp).");
    if (nl[0]!=1 || nl[1]!=*Kp) mexErrMsgTxt("Incompatible dimensions (lkp).");
    lkp0 = (unsigned long long int *)mxGetPr(prhs[5]);


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
        size_t ds[5], i;
        unsigned long long *ptr;
        if (!mxIsNumeric(prhs[9]) || mxIsComplex(prhs[9]) ||
              mxIsSparse(prhs[9]) || !mxIsUint64(prhs[9]))
            mexErrMsgTxt("Skip data must be numeric, real, full and UInt64.");
        nd  = copy_dims(prhs[9],ds);
        if (nd>2) mexErrMsgTxt("Wrong number of dimensions (skip).");
        if (ds[0]*ds[1] > 3) mexErrMsgTxt("Wrong number of elements (skip).");
        ptr = (unsigned long long *)mxGetPr(prhs[9]);
        for(i=0; i<ds[0]*ds[1]; i++)
        {
            if(ptr[i] > 0) skip[i] = ptr[i];
        }
    }

    /* lkp - allocated memory needs freeing later on */
    *lkpp = (size_t *)mxCalloc(sizeof(size_t),*Kp);
    for(k=0; k<*Kp; k++) (*lkpp)[k] = lkp0[k]-1;
}

static void resp_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    static mwSize scalar_dim[2] = {1,1};
    mwSize K, k, nf[5], nm[5], nr[5], skip[3];
    double *m, *b, *W, *nu, *gam, *ll;
    float *mu, *mf, *vf, *r;
    size_t *lkp;

    if (nrhs<9 || nrhs>10 || nlhs>2) mexErrMsgTxt("Incorrect usage");

    parse_rhs(nrhs, prhs, &K, &m, &b, &W, &nu, &gam, &lkp, nm, &mu, nf, &mf, &vf, skip);
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

    ll[0] = call_responsibilities(nf,skip,mf,vf, K,m,b,W,nu,gam, nm[3],lkp,mu, r);
    mxFree(lkp);
}

static void moments_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t K, K1, k, nf[5], nm[5], skip[3], *lkp, dm0, dm1, dm2;
    double *m, *b, *W, *nu, *gam, *s0, *s1, *s2, *ll;
    float *mu, *mf, *vf, *r;

    if (nrhs<9 || nrhs>10 || nlhs>4) mexErrMsgTxt("Incorrect usage");

    parse_rhs(nrhs, prhs, &K, &m, &b, &W, &nu, &gam, &lkp, nm, &mu, nf, &mf, &vf, skip);

    space_needed(nf[3], K, &dm0, &dm1, &dm2);
    plhs[0] = mxCreateDoubleMatrix(dm0,1, mxREAL); s0 = (double *)mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(dm1,1, mxREAL); s1 = (double *)mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(dm2,1, mxREAL); s2 = (double *)mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(  1,1, mxREAL); ll = (double *)mxGetPr(plhs[3]);

    ll[0]   = call_suffstats_missing(nf,mf,vf, K,m,b,W,nu,gam, nm,skip,lkp,mu, s0,s1,s2);

    mxFree(lkp);
}

static void inugrads_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize nf[5];
    size_t K, K1, k, nm[5], dc[5], skip[3], *lkp, dm0, dm1, dm2, c, nd;
    double *m, *b, *W, *nu, *gam, *ll;
    float *mu, *mf, *vf, *g1, *g2;

    if (nrhs!=11 || nlhs>3) mexErrMsgTxt("Incorrect usage");

    parse_rhs(nrhs, prhs, &K, &m, &b, &W, &nu, &gam, &lkp, nm, &mu, nf, &mf, &vf, skip);

    /* c */
    if (!mxIsNumeric(prhs[10]) || mxIsComplex(prhs[10]) ||
          mxIsSparse(prhs[10]) || !mxIsUint64(prhs[10]))
        mexErrMsgTxt("Index must be numeric, real, full and UInt64.");
    nd  = copy_dims(prhs[10],dc);
    if (nd>2 || dc[0]!=1 || dc[1]!=1) mexErrMsgTxt("Index not a scalar.");
    c = ((size_t *)mxGetPr(prhs[10]))[0] - 1;
    if (c<0 || c>=nf[3]) mexErrMsgTxt("Index out of range.");

    plhs[0] = mxCreateNumericArray(3,nf, mxSINGLE_CLASS, mxREAL); g1 = (float *)mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(3,nf, mxSINGLE_CLASS, mxREAL); g2 = (float *)mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    ll      = (double *)mxGetPr(plhs[2]);
    ll[0]   = call_INUgrads(nf,mf,vf, K,m,b,W,nu,gam, nm,skip,lkp,mu, c, g1,g2);
    mxFree(lkp);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    spm_set_num_threads(spm_get_num_threads());

    if ((nrhs>=1) && mxIsChar(prhs[0]))
    {
        int buflen;
        char *fnc_str;
        buflen = mxGetNumberOfElements(prhs[0]);
        fnc_str = (char *)mxCalloc(buflen+1,sizeof(mxChar));
        mxGetString(prhs[0],fnc_str,buflen+1);

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

