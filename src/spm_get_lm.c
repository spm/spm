/*
 * $Id: spm_get_lm.c 6534 2015-08-24 16:02:56Z guillaume $
 * Jesper Andersson
 */

/****************************************************************
 **
 ** Routine that identifies which voxels in a list of coordinates
 ** that are local maxima, and returns a list of indices into
 ** the coordinate list for those maxima.
 **
 ***************************************************************/

#include <math.h>
#include <limits.h>
#include <string.h>
#include "mex.h"

/* Macros */

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif


/* get_index */

mwSignedIndex get_index(mwIndex x, mwIndex y, mwIndex z, mwSize dim[3])
{
   if (x < 1 || x > dim[0] || y < 1 || y > dim[1] || z < 1 || z > dim[2])
   {
       return(-1);
   }
   else
   {
       return((z-1)*dim[0]*dim[1]+(y-1)*dim[0]+x-1);
   }
}

/* is_maxima */

int is_maxima(double        *v,
              mwSize        dim[3],
              mwIndex       x,
              mwIndex       y,
              mwIndex       z,
              unsigned int  cc)
{
   mwSignedIndex ii = 0, i = 0;
   double  cv = 0.0;

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


/* get_maxima */

unsigned int get_maxima(double        *vol,
                        mwSize        vdim[3],
                        double        *list,
                        mwSize        nlist,
                        unsigned int  cc,
                        mwIndex       **ldx)
{
   mwIndex  i = 0, j = 0;
   mwIndex  ix = 0, iy = 0, iz = 0;
   mwSize   ldx_sz = 0;
   mwIndex  ldx_n = 0;

   *ldx = (mwIndex *) mxCalloc((ldx_sz = 1024),sizeof(mwIndex));

   for (i=0, j=0; i<nlist; i++, j+=3)
   {
      /* 
      ** Casting of double to int isn't properly defined in C
      ** (i.e. wether it results in truncation or rounding), 
      ** hence I add a small offset (0.1) to make sure it
      ** works either way.
      */

      ix = ((mwIndex) (list[j]+0.1));
      iy = ((mwIndex) (list[j+1]+0.1));
      iz = ((mwIndex) (list[j+2]+0.1));
      
      if (get_index(ix,iy,iz,vdim) >= 0)
      {
         if (is_maxima(vol,vdim,ix,iy,iz,cc))
         {
            if (ldx_n >= ldx_sz)
            {
               *ldx = (mwIndex *) mxRealloc(*ldx,(ldx_sz += 1024)*sizeof(mwIndex));
            }
            (*ldx)[ldx_n] = i+1;
            ldx_n++;   
         }
      }
   }
   return(ldx_n); 
}


/* Gateway function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   mwIndex       i = 0, j = 0, k = 0, m = 0, n = 0;
   int           tmpint = 0;
   const mwSize  *pdim = NULL;
   mwSize        ndim = 0;
   mwSize        vdim[3];
   mwSize        ln = 0, lm = 0;
   mwSize        n_lindex = 0;
   mwIndex       *lindex = NULL;
   unsigned int  conn;
   double        *vol = NULL;
   double        *lp = NULL;
   double        *list = NULL;
   double        *plindex = NULL;
   
   if (nrhs < 1) mexErrMsgTxt("Not enough input arguments.");
   if (nrhs > 3) mexErrMsgTxt("Too many input arguments.");
   if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");

   /* Get binary map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("Input 'vol' must be numeric, real, full and double.");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if (ndim != 3 && ndim != 2)
   {
      mexErrMsgTxt("Input 'vol' must 2- or 3-dimensional.");
   }
   pdim = mxGetDimensions(prhs[0]);
   vdim[0] = pdim[0]; vdim[1] = pdim[1]; 
   if (ndim == 2) { vdim[2] = 1; } else { vdim[2] = pdim[2]; } 
   vol = mxGetPr(prhs[0]);

   /* Get list of coordinates */
   
   if (nrhs < 2)
   {
      lm = 3; ln = vdim[0] * vdim[1] * vdim[2];
      list = (double *) mxCalloc(lm*ln,sizeof(double));
      for (k=0,m=0,n=0; k<vdim[2]; k++)
      {
         for (j=0; j<vdim[1]; j++)
         {
            for (i=0; i<vdim[0]; i++)
            {
               /* ignore NaNs */
               if (!mxIsNaN(vol[n++]))
               {
                  list[m++] = (double)i + 1.0;
                  list[m++] = (double)j + 1.0;
                  list[m++] = (double)k + 1.0;
               }
               else
               {
                  list[m++] = -1.0;
                  list[m++] = -1.0;
                  list[m++] = -1.0;
               }
            }
         }
      }
   }
   else
   {
      if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
      {
         mexErrMsgTxt("Input 'list' must be numeric, real, full and double.");
      }
      lm = mxGetM(prhs[1]);
      ln = mxGetN(prhs[1]);
      if (!((lm==3 && ndim==3) || (lm==3 && ndim==2) || (lm==2 && ndim==2)))
      {
         mexErrMsgTxt("Input 'list' must be 3xn (or 2xn) list of voxel coordinates.");
      }
      lp = mxGetPr(prhs[1]);
      if (lm==3 && ndim==2)  /* Make sure all z-coordinates equal 1 if 3xn list used with 2D map. */
      {
         for (i=0, j=0; i<ln; i++, j+=3)
         {
            tmpint = ((int) lp[j+2]+0.1);
            if (tmpint != 1)
            {
               mexErrMsgTxt("Input 'list' z-coordinate must be 1 when using 3xn list with 2D map.");
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
   }

   /* Get connectivity criterion */
   
   if (nrhs < 3)
   {
      conn = 18;
   }
   else
   {
      if (!mxIsNumeric(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1)
      {
         mexErrMsgTxt("Input 'n' must be a scalar.");
      }
      conn = (unsigned int)(mxGetScalar(prhs[2]) + 0.1);
   }
   if ((conn!=6) && (conn!=18) && (conn!=26))
   {
      mexErrMsgTxt("Input 'n' must be 6, 18 or 26.");
   }
   
   /* Find list if indices to local maxima. */

   n_lindex = get_maxima(vol,vdim,list,ln,conn,&lindex);

   /* Turn indices into doubles in a MATLAB array. */

   plhs[0] = mxCreateDoubleMatrix(1,n_lindex,mxREAL);
   plindex = mxGetPr(plhs[0]);
   for (i=0; i<n_lindex; i++) {plindex[i] = ((double) lindex[i]);}

   mxFree(list);
   mxFree(lindex);
}
