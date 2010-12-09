#include "mex.h"
#include <stdlib.h>
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

void restore_ramps(double          *mask,
                   unsigned int    dim[3],
                   double          ramp[3],
                   double          *pm);



void restore_ramps(double          *mask,
                   unsigned int    dim[3],
                   double          ramp[3],
                   double          *pm)
{
   int    i=0, j=0, k=0;
   int    ii=0;

   for (i=0; i<dim[0]; i++)
   {
      for (j=0; j<dim[1]; j++)
      {
         for (k=0; k<dim[2]; k++)
     {
        if (mask[ii=index(i,j,k,dim)])
        {
               pm[ii] += ramp[0] * ((double) (i - ((double) (dim[0]-1.0))/2.0));
           pm[ii] += ramp[1] * ((double) (j - ((double) (dim[1]-1.0))/2.0));
           pm[ii] += ramp[2] * ((double) (k - ((double) (dim[2]-1.0))/2.0));
            }
         }
      }
   }

   return;
}


/* Gateway function with error check. */

void mexFunction(int             nlhs,      /* No. of output arguments */
                 mxArray         *plhs[],   /* Output arguments. */ 
                 int             nrhs,      /* No. of input arguments. */
                 const mxArray   *prhs[])   /* Input arguments. */
{
   int            ndim, mask_ndim;
   int            n, i;
   int            ramp_m=0, ramp_n=0;
   const int      *cdim = NULL, *mask_cdim = NULL;
   unsigned int   dim[3];
   double         *mask = NULL;
   double         *pm = NULL;
   double         *opm = NULL;
   double         *ramps = NULL;


   if (nrhs == 0) mexErrMsgTxt("usage: pm = pm_restore_ramp(pm,mask,ramps)");
   if (nrhs != 3) mexErrMsgTxt("pm_restore_ramp: 3 input arguments required");
   if (nlhs != 1) mexErrMsgTxt("pm_restore_ramp: 1 output argument required");

   /* Get phase map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("pm_restore_ramp: pm must be numeric, real, full and double");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if ((ndim < 2) | (ndim > 3))
   {
      mexErrMsgTxt("pm_restore_ramp: pm must be 2 or 3-dimensional");
   }
   cdim = mxGetDimensions(prhs[0]);
   pm = mxGetPr(prhs[0]);

   /* Get mask. */

   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("pm_restore_ramp: mask must be numeric, real, full and double");
   }
   mask_ndim = mxGetNumberOfDimensions(prhs[1]);
   if (mask_ndim != ndim)
   {
      mexErrMsgTxt("pm_restore_ramp: pm and mask must have same dimensionality");
   }
   mask_cdim = mxGetDimensions(prhs[1]);
   for (i=0; i<ndim; i++)
   {
      if (cdim[i] != mask_cdim[i])
      {
         mexErrMsgTxt("pm_estimate_ramp: pm and mask must have same size");
      }
   }
   mask = mxGetPr(prhs[1]);

   /* Fix dimensions to allow for 2D and 3D data. */

   dim[0]=cdim[0]; dim[1]=cdim[1];
   if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 
   for (i=0, n=1; i<ndim; i++)
   {
      n *= dim[i];
   }

   /* Get 3x1 or 1x3 ramp array */

   ramp_m = mxGetM(prhs[2]);
   ramp_n = mxGetN(prhs[2]);
   if (!((ramp_m==1 && ramp_n==3) || (ramp_m==3 && ramp_n==1)))
   {
      mexErrMsgTxt("pm_restore_ramp: ramps must be 3x1 or 1x3");
   }
   ramps = mxGetPr(prhs[2]);

   /* Allocate output phasemap with ramps restored. */

   plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                  mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
   opm = mxGetPr(plhs[0]);

   memcpy(opm,pm,n*sizeof(double));
   restore_ramps(mask,dim,ramps,opm);   
   
   return;
}







