#include "mex.h"
#include <math.h>
#include <limits.h>
#include <string.h>

/* Silly little macros. */

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

void ff_unwrap(double        *ipm,
               double        *vm,
               double        *iwm,
               double        *mask,
               mwSize        dim[3],
               double        thres,
               double        *opm,
               double        *owm);

unsigned int ff_unwrap_dilate(double        *ipm,
                              double        *vm,
                              double        *iwm,
                              double        *mask,
                              mwSize        dim[3],
                              double        thres,
                              double        *opm,
                              double        *owm);

double neighbour(mwIndex       i,
                 mwIndex       j,
                 mwIndex       k,
                 mwSize        dim[3],
                 double        *wm,
                 double        *pm,
                 int           cc);

void ff_unwrap(double        *ipm,
               double        *vm,
               double        *iwm,
               double        *mask,
               mwSize        dim[3],
               double        thres,
               double        *opm,
               double        *owm)
{
    double *tpm[2] = {NULL,NULL};
    double *twm[2] = {NULL,NULL};
    double *tmp = NULL;
    
    if (ff_unwrap_dilate(ipm,vm,iwm,mask,dim,thres,opm,owm))
    {
        tpm[0] = opm;
        twm[0] = owm;
        tpm[1] = (double *) mxCalloc(dim[0]*dim[1]*dim[2],sizeof(double));
        twm[1] = (double *) mxCalloc(dim[0]*dim[1]*dim[2],sizeof(double));
        while (ff_unwrap_dilate(tpm[0],vm,twm[0],mask,dim,thres,tpm[1],twm[1]))
        {
            tmp = tpm[0];
            tpm[0] = tpm[1];
            tpm[1] = tmp;
            tmp = twm[0];
            twm[0] = twm[1];
            twm[1] = tmp;
        }
        if (tpm[1]==opm)
        {
            mxFree(tpm[0]);
            mxFree(twm[0]);
        }
        else
        {
            memcpy(opm,tpm[1],dim[0]*dim[1]*dim[2]*sizeof(double));
            memcpy(owm,twm[1],dim[0]*dim[1]*dim[2]*sizeof(double));
            mxFree(tpm[1]);
            mxFree(twm[1]);
        }
    }
    
    return;
}

unsigned int ff_unwrap_dilate(double        *ipm,
                              double        *vm,
                              double        *iwm,
                              double        *mask,
                              mwSize        dim[3],
                              double        thres,
                              double        *opm,
                              double        *owm)
{
    mwIndex       i=0, j=0, k=0;
    int           ii = 0;
    unsigned int  nwrapped = 0;
    double        np = 0;
    
    for (i=0; i<dim[0]; i++)
    {
        for (j=0; j<dim[1]; j++)
        {
            for (k=0; k<dim[2]; k++)
            {
                ii=index(i,j,k,dim);
                if (mask[ii] && !iwm[ii] && vm[ii] < thres)
                {
                    if ((np = neighbour(i,j,k,dim,iwm,ipm,6)))
                    {
                        owm[ii] = 1.0;
                        opm[ii] = ipm[ii];
                        while ((opm[ii] - np) > PI) {opm[ii] -= 2*PI;}
                        while ((opm[ii] - np) < -PI) {opm[ii] += 2*PI;}
                        nwrapped++;
                    }
                    else
                    {
                        owm[ii] = iwm[ii];
                        opm[ii] = ipm[ii];
                    }
                }
                else
                {
                    owm[ii] = iwm[ii];
                    opm[ii] = ipm[ii];
                }
            }
        }
    }
    
    return(nwrapped);
}

double neighbour(mwIndex       i,
                 mwIndex       j,
                 mwIndex       k,
                 mwSize        dim[3],
                 double        *wm,
                 double        *pm,
                 int           cc)
{
    double    np = 0.0;
    int       nn = 0;
    int       ii = 0;
    
    if (cc >= 6)
    {
        if (i && wm[ii=index(i-1,j,k,dim)]) {np+=pm[ii]; nn++;}
        if (j && wm[ii=index(i,j-1,k,dim)]) {np+=pm[ii]; nn++;}
        if (k && wm[ii=index(i,j,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (i<(dim[0]-1) && wm[ii=index(i+1,j,k,dim)]) {np+=pm[ii]; nn++;}
        if (j<(dim[1]-1) && wm[ii=index(i,j+1,k,dim)]) {np+=pm[ii]; nn++;}
        if (k<(dim[2]-1) && wm[ii=index(i,j,k+1,dim)]) {np+=pm[ii]; nn++;}
    }
    if (cc >= 18)
    {
        if (i && j && wm[ii=index(i-1,j-1,k,dim)]) {np+=pm[ii]; nn++;}
        if (i && k && wm[ii=index(i-1,j,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (j && k && wm[ii=index(i,j-1,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (i && j<(dim[1]-1) && wm[ii=index(i-1,j+1,k,dim)]) {np+=pm[ii]; nn++;}
        if (j && i<(dim[0]-1) && wm[ii=index(i+1,j-1,k,dim)]) {np+=pm[ii]; nn++;}
        if (i && k<(dim[2]-1) && wm[ii=index(i-1,j,k+1,dim)]) {np+=pm[ii]; nn++;}
        if (k && i<(dim[0]-1) && wm[ii=index(i+1,j,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (j && k<(dim[2]-1) && wm[ii=index(i,j-1,k+1,dim)]) {np+=pm[ii]; nn++;}
        if (k && j<(dim[1]-1) && wm[ii=index(i,j+1,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (i<(dim[0]-1) && j<(dim[1]-1) && wm[ii=index(i+1,j+1,k,dim)]) {np+=pm[ii]; nn++;}
        if (i<(dim[0]-1) && k<(dim[2]-1) && wm[ii=index(i+1,j,k+1,dim)]) {np+=pm[ii]; nn++;}
        if (j<(dim[1]-1) && k<(dim[2]-1) && wm[ii=index(i,j+1,k+1,dim)]) {np+=pm[ii]; nn++;}
    }
    if (cc >= 26)
    {
        if (i && j && k && wm[ii=index(i-1,j-1,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (i && j && k<(dim[2]-1) && wm[ii=index(i-1,j-1,k+1,dim)]) {np+=pm[ii]; nn++;}
        if (i && j<(dim[1]-1) && k && wm[ii=index(i-1,j+1,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (i<(dim[0]-1) && j && k && wm[ii=index(i+1,j-1,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (i && j<(dim[1]-1) && k<(dim[2]-1) && wm[ii=index(i-1,j+1,k+1,dim)]) {np+=pm[ii]; nn++;}
        if (i<(dim[0]-1) && j && k<(dim[2]-1) && wm[ii=index(i+1,j-1,k+1,dim)]) {np+=pm[ii]; nn++;}
        if (i<(dim[0]-1) && j<(dim[1]-1) && k && wm[ii=index(i+1,j+1,k-1,dim)]) {np+=pm[ii]; nn++;}
        if (i<(dim[0]-1) && j<(dim[1]-1) && k<(dim[2]-1) && wm[ii=index(i+1,j+1,k+1,dim)]) {np+=pm[ii]; nn++;}
    }

   if (nn) return(np/((double) nn));
   else return(0.0);
}

/* Gateway function with error check. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   mwSize         ndim = 0, chk_ndim = 0;
   mwIndex        n = 0, i = 0;
   mwIndex        tm = 0, tn = 0;
   const mwSize   *cdim = NULL, *chk_cdim = NULL;
   mwSize         dim[3];
   double         *pm = NULL;
   double         *vm = NULL;
   double         *wm = NULL;
   double         *mask = NULL;
   double         *thres = NULL;
   double         *opm = NULL;
   double         *owm = NULL;
   double         *tpm = NULL;
   double         *twm = NULL;

   if (nrhs == 0) mexErrMsgTxt("usage: [pm,wm]=pm_ff_unwrap(pn,vm,wm,mask,thres)");
   if (nrhs != 5) mexErrMsgTxt("pm_ff_unwrap: 5 input arguments required");
   if (nlhs != 2) mexErrMsgTxt("pm_ff_unwrap: 2 output argument required");

   /* Get phase map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("pm_ff_unwrap: pm must be numeric, real, full and double");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if ((ndim < 2) | (ndim > 3))
   {
      mexErrMsgTxt("pm_ff_unwrap: pm must be 2 or 3-dimensional");
   }
   cdim = mxGetDimensions(prhs[0]);
   pm = mxGetPr(prhs[0]);

   /* Get variance-map. */

   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("pm_ff_unwrap: vm must be numeric, real, full and double");
   }
   chk_ndim = mxGetNumberOfDimensions(prhs[1]);
   if (chk_ndim != ndim)
   {
      mexErrMsgTxt("pm_ff_unwrap: pm and vm must have same dimensionality");
   }
   chk_cdim = mxGetDimensions(prhs[1]);
   for (i=0; i<ndim; i++)
   {
      if (cdim[i] != chk_cdim[i])
      {
         mexErrMsgTxt("pm_ff_unwrap: pm and vm must have same size");
      }
   }
   vm = mxGetPr(prhs[1]);

   /* Get wrap-map. */

   if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
   {
      mexErrMsgTxt("pm_ff_unwrap: wm must be numeric, real, full and double");
   }
   chk_ndim = mxGetNumberOfDimensions(prhs[2]);
   if (chk_ndim != ndim)
   {
      mexErrMsgTxt("pm_ff_unwrap: pm and wm must have same dimensionality");
   }
   chk_cdim = mxGetDimensions(prhs[2]);
   for (i=0; i<ndim; i++)
   {
      if (cdim[i] != chk_cdim[i])
      {
         mexErrMsgTxt("pm_ff_unwrap: pm and wm must have same size");
      }
   }
   wm = mxGetPr(prhs[2]);

   /* Get mask. */

   if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
   {
      mexErrMsgTxt("pm_ff_unwrap: mask must be numeric, real, full and double");
   }
   chk_ndim = mxGetNumberOfDimensions(prhs[3]);
   if (chk_ndim != ndim)
   {
      mexErrMsgTxt("pm_ff_unwrap: pm and mask must have same dimensionality");
   }
   chk_cdim = mxGetDimensions(prhs[3]);
   for (i=0; i<ndim; i++)
   {
      if (cdim[i] != chk_cdim[i])
      {
         mexErrMsgTxt("pm_ff_unwrap: pm and mask must have same size");
      }
   }
   mask = mxGetPr(prhs[3]);

   /* Fix dimensions to allow for 2D and 3D data. */

   dim[0]=cdim[0]; dim[1]=cdim[1];
   if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 
   for (i=0, n=1; i<ndim; i++)
   {
      n *= dim[i];
   }

   /* Get scalar or array of thresholds. */

   if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) || mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4]))
   {
      mexErrMsgTxt("pm_ff_unwrap: thres must be numeric, real, full and double");
   }
   tm = mxGetM(prhs[4]);
   tn = mxGetN(prhs[4]);
   if (tm>1 && tn>1)
   {
      mexErrMsgTxt("pm_ff_unwrap: thres must be scalar or 1D array");
   }
   tn = MAX(tm,tn);
   thres = mxGetPr(prhs[4]);
      
   /* Allocate memory for output. */

   plhs[0] = mxCreateNumericArray(ndim,dim,mxDOUBLE_CLASS,mxREAL);
   opm = mxGetPr(plhs[0]);
   plhs[1] = mxCreateNumericArray(ndim,dim,mxDOUBLE_CLASS,mxREAL);
   owm = mxGetPr(plhs[1]);

   /* Allocate scratch pad. */
   
   tpm = (double *) mxCalloc(n,sizeof(double));
   twm = (double *) mxCalloc(n,sizeof(double));
   memcpy(tpm,pm,n*sizeof(double));
   memcpy(twm,wm,n*sizeof(double));


   for (i=0; i<tn; i++)
   {
      ff_unwrap(tpm,vm,twm,mask,dim,thres[i],opm,owm);
      memcpy(tpm,opm,n*sizeof(double));
      memcpy(twm,owm,n*sizeof(double));     
   }

   mxFree(tpm);
   mxFree(twm);

   return;
}
