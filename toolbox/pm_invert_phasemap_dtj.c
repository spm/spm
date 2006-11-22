#include "mex.h"
#include <math.h>

#ifndef MAX
#define MAX(a,b)     (((a)>(b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)     (((a)<(b)) ? (a) : (b))
#endif

#ifndef INDX
#define INDX(f,p,s,dim) (((s)-1)*dim[0]*dim[1]+((p)-1)*dim[0]+((f)-1)) 
#endif
 
/* Function doing the job. */

void invert_phasemap_dtj(unsigned int    *dim,
                         double          *pm,
                         double          *ipm)
{
   unsigned int     i = 0;
   unsigned int     oi = 0;
   unsigned int     f = 0;     /* Index in frequency encode direction. */
   unsigned int     p = 0;     /* Index in phase encode direction. */
   unsigned int     s = 0;     /* Index in slice selection direction. */

   for (s=1; s<=dim[2]; s++)    /* New slice. */
   {
      for (f=1; f<=dim[0]; f++)  /* New column in phase encode direction. */
      {
	 oi = 1;
         for (p=1; p<=dim[1]; p++)
	 {
	    for (i=oi; i<=dim[1] && (pm[INDX(f,i,s,dim)]+i)<p; i++) ; /* N.B. */
            if (i>1 && i<=dim[1])
	    {
	       ipm[INDX(f,p,s,dim)] = ((double) i) - 1.0 - ((double) p) + 
                                      (((double) p)-pm[INDX(f,i-1,s,dim)]-((double) (i-1))) / 
                                      (pm[INDX(f,i,s,dim)]-pm[INDX(f,i-1,s,dim)]+1.0);
            }
            else
            {
	       ipm[INDX(f,p,s,dim)] = DBL_MAX;
            }
            oi = MAX(1,i-1);
         }

         for (i=1; i<=dim[1] && ipm[INDX(f,i,s,dim)]==DBL_MAX; i++) ; /* N.B. */
         for (p=i-1; p>0; p--)
	 {
	    ipm[INDX(f,p,s,dim)] = ipm[INDX(f,i,s,dim)];
	 }

         for (i=dim[1]; i>0 && ipm[INDX(f,i,s,dim)]==DBL_MAX; i--) ; /* N.B. */
         for (p=i+1; p<=dim[1]; p++)
	 {
	    ipm[INDX(f,p,s,dim)] = ipm[INDX(f,i,s,dim)];
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
   int            i = 0;
   unsigned int   dim[3] = {0, 0, 0};
   unsigned int   dim_n = 0;
   unsigned int   dim_m = 0;
   unsigned int   pm_ndim = 0;
   unsigned int   pm_dim[3] = {0, 0, 0};
   double         ddim[3] = {0.0, 0.0, 0.0};

   if (nrhs == 0) mexErrMsgTxt("usage: ipm = invert_phasemap_dtj(pm,dim)");
   if (nrhs != 2) mexErrMsgTxt("make_A: 2 input argument required");
   if (nlhs != 1) mexErrMsgTxt("make_A: 1 output arguments required");

   /* Get phasemap (deformation field). */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("invert_phasemap_dtj: pm must be numeric, real, full and double");
   }

   if ((pm_ndim = mxGetNumberOfDimensions(prhs[0])) > 3)
   {
      mexErrMsgTxt("invert_phasemap_dtj: pm must be 2 or 3 dimensional");
   }
   memcpy(pm_dim,mxGetDimensions(prhs[0]),pm_ndim*sizeof(int));
 
   /* Note that dim is allways 3 values   */ 
   /* regardless of dimensionality of pm. */
  
   dim_m = mxGetM(prhs[1]);
   dim_n = mxGetN(prhs[1]);
   if (!((dim_m==1 && dim_n==3) || (dim_m==3 && dim_n==1)))
   {
      mexErrMsgTxt("invert_phsemap_dtj: dim should be a 1x3 or 3x1 vector");
   }
   memcpy(ddim,mxGetPr(prhs[1]),3*sizeof(double));
   for (i=0; i<3; i++) { dim[i] = ((unsigned int) ddim[i]);}

   /* Allocate memory for output. */

   plhs[0] = mxCreateNumericArray(pm_ndim,pm_dim,mxDOUBLE_CLASS,mxREAL);

   invert_phasemap_dtj(dim,mxGetPr(prhs[0]),mxGetPr(plhs[0]));
      
   return;
}





