#ifndef lint
static char sccsid[]="%W% (c) Jesper Andersson %E%";
#endif

#include "mex.h"
#include <math.h>

#ifndef max
#define max(a,b) ((a)>(b) ? (a) : (b))
#endif

/* Function prototypes. */

int connected(unsigned int    i,
              unsigned int    j,
              unsigned int    k,
              unsigned int    dim[3],
              double          *visit);

int index(unsigned int   i,
          unsigned int   j,
          unsigned int   k,
          unsigned int   dim[3]);


/* Function doing the job. */

void uw_dilate(double        *inmap,
               unsigned int  dim[3],
               double        *omap)
{
   unsigned int   i, j, k;
   int            indx;
   
   for (k=1; k<=dim[2]; k++)
   {
      for (j=1; j<=dim[1]; j++)
      {
         for (i=1; i<=dim[0]; i++)
	 {
	   if (!inmap[(indx = index(i,j,k,dim))])
	   {
	      if (connected(i,j,k,dim,inmap) > 0)
	      {
	         omap[indx] = 1;
              }
           }
	 }
      }
   }

   for (i=0; i<dim[0]*dim[1]*dim[2]; i++)
   {
      omap[i] += inmap[i];
   }
               
   return;
}

/* Utility function that returns the index */
/* of a super-threshold voxel connected to */
/* that given by i,j, and k. If there are  */
/* no adjacent super-threshold values it   */
/* returns -1.                             */

int connected(unsigned int    i,
              unsigned int    j,
              unsigned int    k,
              unsigned int    dim[3],
              double          *visit)
{
  int   indx = 0;

  if ((indx = index(i+1,j,k,dim)) > 0 & visit[indx] > 0) return(indx);
  if ((indx = index(i-1,j,k,dim)) > 0 & visit[indx] > 0) return(indx);
  if ((indx = index(i,j+1,k,dim)) > 0 & visit[indx] > 0) return(indx);
  if ((indx = index(i,j-1,k,dim)) > 0 & visit[indx] > 0) return(indx);
  if ((indx = index(i,j,k+1,dim)) > 0 & visit[indx] > 0) return(indx);
  if ((indx = index(i,j,k-1,dim)) > 0 & visit[indx] > 0) return(indx);

  return(-1);
}

/* Utility function that returns index into */
/* 1D array with range checking.            */
 
int index(unsigned int   i,
          unsigned int   j,
          unsigned int   k,
          unsigned int   dim[3])
{
   int   indx;

   if (i<=1 | i>=dim[0] | j<=1 | j>=dim[1] | k<=1 | k>=dim[2]) return(-1);
   else return((k-1)*dim[0]*dim[1]+(j-1)*dim[0]+i-1);
} 


void mexFunction(int             nlhs,      /* No. of output arguments */
                 mxArray         *plhs[],   /* Output arguments. */ 
                 int             nrhs,      /* No. of input arguments. */
                 const mxArray   *prhs[])   /* Input arguments. */
{
   int            i;
   unsigned int   n;
   unsigned int   ndim;
   const int      *indim = NULL;
   unsigned int   dim[3];
   double         *dindim = NULL;
   double         *inmap = NULL;
   double         *outmap = NULL;

   if (nrhs == 0) mexErrMsgTxt("usage: omap=uw_dilate(inmap,dim)");
   if (nrhs<1 | nrhs>2) mexErrMsgTxt("uw_dilate: 1 or 2 input arguments required");
   if (nlhs != 1) mexErrMsgTxt("uw_dilate: 1 output argument required");

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("uw_dilate: inmap must be numeric, real, full and double");
   }

   if (nrhs==1)
   {
      ndim = mxGetNumberOfDimensions(prhs[0]);
      indim = mxGetDimensions(prhs[0]);
      n = mxGetNumberOfElements(prhs[0]);
   }
   else
   {
      if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
      {
         mexErrMsgTxt("uw_dilate: dim must be numeric, real, full and double");
      }
      ndim = max(mxGetM(prhs[1]),mxGetN(prhs[1]));
      dindim = mxGetPr(prhs[1]);
      for (i=0, n=1; i<ndim; i++) n *= ((int) dindim[i]);
      if (n != mxGetNumberOfElements(prhs[0]))
      { 
         mexErrMsgTxt("uw_dilate: size mismatch between inmap and dim");
      }
   }
   inmap = mxGetPr(prhs[0]);

   if (ndim<2 | ndim>3)
   {
      mexErrMsgTxt("uw_dilate: inmap must be 2D or 3D");
   }
   else if (ndim == 2)
   {  
      if (nrhs == 1) { dim[0] = indim[0]; dim[1] = indim[1]; }
      else {dim[0] = ((int) dindim[0]); dim[1] = ((int) dindim[1]); }
      dim[2] = 1;
      ndim == 3;
   }
   else if (ndim == 3)   
   {
      if (nrhs == 1) {dim[0] = indim[0]; dim[1] = indim[1]; dim[2] = indim[2]; }
      else {dim[0] = ((int) dindim[0]); dim[1] = ((int) dindim[1]); dim[2] = ((int) dindim[2]); }
   }

   /* Allocate memory for output. */
   
   plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                  mxGetDimensions(prhs[0]),
                                  mxGetClassID(prhs[0]),mxREAL);
   outmap = mxGetPr(plhs[0]);
   
   /* Initialise output maps to zeros. */

   memset(outmap,0,n*sizeof(double));

   uw_dilate(inmap,dim,outmap);
}

