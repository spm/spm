/*
 * $Id: spm_get_lm.c 1790 2008-06-05 11:27:02Z spm $
 * Jesper Andersson
 */

/****************************************************************
 **
 ** Routine that identifies which voxels in a list of coordinates
 ** that are local maxima, and returns a list of indicies into
 ** the coordinate list for those maxima.
 **
 ***************************************************************/

#include <math.h>
#include <limits.h>
#include <string.h>
#include "mex.h"

/* Silly little macros. */

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif


/* Function prototypes. */

unsigned int get_maxima(double        *vol,
                        unsigned int  vdim[3],
                        double        *list,
                        unsigned int  nlist,
                        unsigned int  cc,
                        unsigned int  **lindex);

int get_index(int       x,
          int           y,
          int           z,
          unsigned int  dim[3]);

int is_maxima(double        *v,
              unsigned int  vdim[3],
              int           x,
              int           y,
              int           z,
              unsigned int  cc);

/* Here starts actual code. */

unsigned int get_maxima(double        *vol,
                        unsigned int  vdim[3],
                        double        *list,
                        unsigned int  nlist,
                        unsigned int  cc,
                        unsigned int  **ldx)
{
   int           i = 0, j = 0;
   int           ix = 0, iy = 0, iz = 0;
   unsigned int  ldx_sz = 0, ldx_n = 0;

   *ldx = (unsigned int *) mxCalloc((ldx_sz = 1000),sizeof(unsigned int));

   for (i=0, j=0; i<nlist; i++, j+=3)
   {
      /* 
      ** Casting of double to int isn't properly defined in C
      ** (i.e. wether it results in truncation or rounding), 
      ** hence I add a small offset (0.1) to make sure it
      ** works either way.
      */

      ix = ((int) (list[j]+0.1)); iy = ((int) (list[j+1]+0.1)); iz = ((int) (list[j+2]+0.1));
      
      if (get_index(ix,iy,iz,vdim) > 0)
      {
         if (is_maxima(vol,vdim,ix,iy,iz,cc))
         {
            if (ldx_n >= ldx_sz)
            {
               *ldx = (unsigned int *) mxRealloc(*ldx,(ldx_sz += 1000)*sizeof(unsigned int));
            }
            (*ldx)[ldx_n] = i+1;
            ldx_n++;   
         }
      }
   }
   return(ldx_n); 
}


int get_index(int       x,
          int           y,
          int           z,
          unsigned int  dim[3])
{
   if (x < 1 || x > dim[0] || y < 1 || y > dim[1] || z < 1 || z > dim[2]) {return(-1);}
   else {return((z-1)*dim[0]*dim[1]+(y-1)*dim[0]+x-1);}
}


int is_maxima(double        *v,
              unsigned int  dim[3],
              int           x,
              int           y,
              int           z,
              unsigned int  cc)
{
   int   ii = 0, i = 0;
   double cv = 0.0;

   if ((ii=get_index(x,y,z,dim))<0) {return(0);}
   cv = v[ii];

   if (cc >= 6)
   {
      if ((i=get_index(x+1,y,z,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y,z,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x,y+1,z,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x,y-1,z,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x,y,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x,y,z-1,dim))>0 && v[i] > cv) {return(0);}
   }
   if (cc >= 18)
   {
      if ((i=get_index(x+1,y+1,z,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x+1,y-1,z,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x+1,y,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x+1,y,z-1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y+1,z,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y-1,z,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y,z-1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x,y+1,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x,y+1,z-1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x,y-1,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x,y-1,z-1,dim))>0 && v[i] > cv) {return(0);}
   }
   if (cc == 26)
   {
      if ((i=get_index(x+1,y+1,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x+1,y+1,z-1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x+1,y-1,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x+1,y-1,z-1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y+1,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y+1,z-1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y-1,z+1,dim))>0 && v[i] > cv) {return(0);}
      if ((i=get_index(x-1,y-1,z-1,dim))>0 && v[i] > cv) {return(0);}
   }
     
   return(1);
}

/* Gateway function with error check. */

void mexFunction(int            nlhs,      /* No. of output arguments */
                 mxArray        *plhs[],   /* Output arguments. */ 
                 int            nrhs,      /* No. of input arguments. */
                 const mxArray  *prhs[])   /* Input arguments. */
{
   int                 i = 0, j = 0, k = 0;
   int                 tmpint = 0;
   const int           *pdim = NULL;
   unsigned int        ndim = 0;
   unsigned int        vdim[3];
   unsigned int        ln = 0, lm = 0;
   unsigned int        n_lindex = 0;
   unsigned int        *lindex = NULL;
   double              *vol = NULL;
   double              *lp = NULL;
   double              *list = NULL;
   double              *plindex = NULL;
   
   if (nrhs < 2) mexErrMsgTxt("Not enough input arguments.");
   if (nrhs > 2) mexErrMsgTxt("Too many input arguments.");
   if (nlhs < 1) mexErrMsgTxt("Not enough output arguments");
   if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");

   /* Get binary map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("spm_get_lm: VOL must be numeric, real, full and double");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if (ndim != 3 && ndim != 2)
   {
      mexErrMsgTxt("spm_get_lm: VOL must 2- or 3-dimensional");
   }
   pdim = mxGetDimensions(prhs[0]);
   vdim[0] = pdim[0]; vdim[1] = pdim[1]; 
   if (ndim == 2) {vdim[2]=1;} else {vdim[2]=pdim[2];} 
   vol = mxGetPr(prhs[0]);

   /* Get list of coordinates */
   
   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("spm_get_lm: L must be numeric, real, full and double");
   }
   lm = mxGetM(prhs[1]);
   ln = mxGetN(prhs[1]);
   if (!((lm==3 && ndim==3) || (lm==3 && ndim==2) || (lm==2 && ndim==2)))
   {
      mexErrMsgTxt("spm_get_lm: L must be 3xn (or 2xn) list of voxel coordinates");
   }
   lp = mxGetPr(prhs[1]);
   if (lm==3 && ndim==2)  /* Make sure all z-coordinates equal 1 if 3xn list used with 2D map. */
   {
      for (i=0, j=0; i<ln; i++, j+=3)
      {
         tmpint = ((int) lp[j+2]+0.1);
         if (tmpint != 1)
         {
            mexErrMsgTxt("spm_get_lm: z-coordinate must be 1 when using 3xn list with 2D map");
         }
      }
   }
   list = (double *) mxCalloc(3*ln,sizeof(double));
   if (lm==2) /* Extend to 3xn list if 2xn list was supplied with 2D-map. */
   {
      for (i=0, j=0, k=0; i<ln; i++, j+=3, k+=2)
      {
         list[j] = lp[k]; list[j+1] = lp[k+1]; list[j+2] = 1.0; 
      }
   }
   else
   {
      memcpy(list,lp,3*ln*sizeof(double));
   }

   /* Find list if indicies to local maxima. */

   n_lindex = get_maxima(vol,vdim,list,ln,18,&lindex);

   /* Turn indicies into doubles in a Matlab array. */

   plhs[0] = mxCreateDoubleMatrix(1,n_lindex,mxREAL);
   plindex = mxGetPr(plhs[0]);
   for (i=0; i<n_lindex; i++) {plindex[i] = ((double) lindex[i]);}

   mxFree(list);
   mxFree(lindex);

   return;
}
