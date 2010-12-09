#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <limits.h>

int getindex(int            i,
          int            j,
          int            k,
          unsigned int   dim[3]);

void smooth(double         *pm,
            double         *wmap,
            unsigned int   dim[3],
            double         *krnl,
            unsigned int   kdim[3],
            double         *opm);


void smooth(double         *pm,
            double         *wmap,
            unsigned int   dim[3],
            double         *krnl,
            unsigned int   kdim[3],
            double         *opm)
{
   int            i=0, j=0, k=0;
   int            ki=0, kj=0, kk=0;
   int            ndx=0, kndx=0;
   double         ii=0.0, wgt=0.0, twgt=0.0;
  
   for (i=0; i<dim[0]; i++)
   {
      for (j=0; j<dim[1]; j++)
      {
     for (k=0; k<dim[2]; k++)
     {
        ndx = getindex(i,j,k,dim);
        twgt = 0.0;
            ii = 0.0;
        for (ki=0; ki<kdim[0]; ki++)
        {
           for (kj=0; kj<kdim[1]; kj++)
           {
              for (kk=0; kk<kdim[2]; kk++)
              {
                     kndx = getindex(i-(kdim[0]/2)+ki,j-(kdim[1]/2)+kj,k-(kdim[2]/2)+kk,dim);
             if (kndx > -1)
             {
                wgt = krnl[getindex(ki,kj,kk,kdim)] * wmap[kndx];
                        ii += pm[kndx] * wgt;
                        twgt += wgt;
                     }
                  }
               }
            }
            if (twgt)
        {
           opm[ndx] = ii/twgt;
            }
            else
        {
               opm[ndx]=pm[ndx];
            }
         }
      }
   }
   return;
}

/* Utility function that returns index into */
/* 1D array with range checking.            */
 
int getindex(int            i,
          int            j,
          int            k,
          unsigned int   dim[3])
{
   if (i<0 | i>(dim[0]-1) | j<0 | j>(dim[1]-1) | k<0 | k>(dim[2]-1)) return(-1);
   else return(k*dim[0]*dim[1]+j*dim[0]+i);
}


/* Gateway function with error check. */

void mexFunction(int             nlhs,      /* No. of output arguments */
                 mxArray         *plhs[],   /* Output arguments. */ 
                 int             nrhs,      /* No. of input arguments. */
                 const mxArray   *prhs[])   /* Input arguments. */
{
   int            ndim, wmap_ndim, krn_ndim;
   int            n, i;
   const int      *cdim = NULL, *wmap_cdim = NULL, *krn_cdim = NULL;
   unsigned int   dim[3], kdim[3];
   double         *pm = NULL;
   double         *wmap = NULL;
   double         *opm = NULL;
   double         *krnl = NULL;


   if (nrhs == 0) mexErrMsgTxt("usage: pm = pm_pad(pm,wmap,kernel)");
   if (nrhs != 3) mexErrMsgTxt("pm_smooth_phasemap_dtj: 3 input arguments required");
   if (nlhs != 1) mexErrMsgTxt("pm_smooth_phasemap_dtj: 1 output argument required");

   /* Get phase map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("pm_smooth_phasemap_dtj: pm must be numeric, real, full and double");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if ((ndim < 2) | (ndim > 3))
   {
      mexErrMsgTxt("pm_smooth_phasemap_dtj: pm must be 2 or 3-dimensional");
   }
   cdim = mxGetDimensions(prhs[0]);
   pm = mxGetPr(prhs[0]);

   /* Get weight-map (reciprocal of variance). */

   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("pm_smooth_phasemap_dtj: wmap must be numeric, real, full and double");
   }
   wmap_ndim = mxGetNumberOfDimensions(prhs[1]);
   if (wmap_ndim != ndim)
   {
      mexErrMsgTxt("pm_smooth_phasemap_dtj: pm and wmap must have same dimensionality");
   }
   wmap_cdim = mxGetDimensions(prhs[1]);
   for (i=0; i<ndim; i++)
   {
      if (cdim[i] != wmap_cdim[i])
      {
         mexErrMsgTxt("pm_smooth_phasemap_dtj: pm and wmap must have same size");
      }
   }
   wmap = mxGetPr(prhs[1]);

   /* Fix dimensions to allow for 2D and 3D data. */

   dim[0]=cdim[0]; dim[1]=cdim[1];
   if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 
   for (i=0, n=1; i<ndim; i++)
   {
      n *= dim[i];
   }

   /* Get kernel */

   if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
   {
      mexErrMsgTxt("pm_smooth_phasemap_dtj: kernel must be numeric, real, full and double");
   }
   krn_ndim = mxGetNumberOfDimensions(prhs[2]);
   if (krn_ndim != ndim)
   {
      mexErrMsgTxt("pm_smooth_phasemap_dtj: pm and kernel must have same dimensionality");
   }
   krn_cdim = mxGetDimensions(prhs[2]);
   krnl = mxGetPr(prhs[2]);
   kdim[0]=krn_cdim[0]; kdim[1]=krn_cdim[1];
   if (krn_ndim==2) {kdim[2]=1; krn_ndim=3;} else {kdim[2]=krn_cdim[2];} 

   /* Allocate mem for smoothed output phasemap. */

   plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                  mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
   opm = mxGetPr(plhs[0]);

   smooth(pm,wmap,dim,krnl,kdim,opm);
   
   return;
}


