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

double wrap(double  angle);

void estimate_ramps(double          *pm,
                    double          *mask,
                    unsigned int    dim[3],
                    double          ramp[3]);

void remove_ramps(double          *mask,
                  unsigned int    dim[3],
                  double          ramp[3],
                  double          *pm);



double wrap(double  angle)
{
  while (angle > PI) {angle -= 2*PI;}
  while (angle < -PI) {angle += 2*PI;}

  return(angle);
}

void estimate_ramps(double          *pm,
                    double          *mask,
                    unsigned int    dim[3],
                    double          ramp[3])
{
   int    i=0, j=0, k=0;
   int    nx=0, ny=0, nz=0;
   int    c=0, ci=0, cj=0, ck=0;

   for (i=0; i<dim[0]; i++)
   {
      for (j=0; j<dim[1]; j++)
      {
         for (k=0; k<dim[2]; k++)
     {
        if (mask[c=index(i,j,k,dim)])
        {
           if ((i<(dim[0]-1)) && mask[ci=index(i+1,j,k,dim)])
           {
          ramp[0] += wrap(pm[ci]-pm[c]);
                  nx++; 
               }
           if ((j<(dim[1]-1)) && mask[cj=index(i,j+1,k,dim)])
           {
          ramp[1] += wrap(pm[cj]-pm[c]);
                  ny++;
               }
           if ((k<(dim[2]-1)) && mask[ck=index(i,j,k+1,dim)])
           {
          ramp[2] += wrap(pm[ck]-pm[c]);
                  nz++;
               }
            }
         }
      }
   }

   ramp[0] /= ((double) nx);
   ramp[1] /= ((double) ny);
   if (nz) ramp[2] /= ((double) nz);
   else ramp[2] = 0.0;

   return;
}

void remove_ramps(double          *mask,
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
               pm[ii] -= ramp[0] * ((double) (i - ((double) (dim[0]-1.0))/2.0));
           pm[ii] -= ramp[1] * ((double) (j - ((double) (dim[1]-1.0))/2.0));
           pm[ii] -= ramp[2] * ((double) (k - ((double) (dim[2]-1.0))/2.0));
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
   const int      *cdim = NULL, *mask_cdim = NULL;
   unsigned int   dim[3];
   double         *mask = NULL;
   double         *pm = NULL;
   double         *opm = NULL;
   double         *ramps = NULL;


   if (nrhs == 0) mexErrMsgTxt("usage: [ramp,pm] = pm_estimate_ramp(pm,mask)");
   if (nrhs != 2) mexErrMsgTxt("pm_estimate_ramp: 2 input arguments required");
   if (nlhs != 2) mexErrMsgTxt("pm_estimate_ramp: 2 output arguments required");

   /* Get phase map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("pm_estimate_ramp: pm must be numeric, real, full and double");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if ((ndim < 2) | (ndim > 3))
   {
      mexErrMsgTxt("pm_estimate_ramp: pm must be 2 or 3-dimensional");
   }
   cdim = mxGetDimensions(prhs[0]);
   pm = mxGetPr(prhs[0]);

   /* Get mask. */

   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("pm_estimate_ramp: pm must be numeric, real, full and double");
   }
   mask_ndim = mxGetNumberOfDimensions(prhs[1]);
   if (mask_ndim != ndim)
   {
      mexErrMsgTxt("pm_estimate_ramp: pm and mask must have same dimensionality");
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

   /* Allocate output 3x1 matrix. */

   plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
   ramps = mxGetPr(plhs[0]);

   /* Allocate output phasemap without ramps. */

   plhs[1] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                  mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
   opm = mxGetPr(plhs[1]);

   estimate_ramps(pm,mask,dim,ramps);
   memcpy(opm,pm,n*sizeof(double));
   remove_ramps(mask,dim,ramps,opm);   
   
   return;
}







