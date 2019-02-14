#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>


#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif


unsigned int make_vectors(double        *rima,
                          double        *pm,
                          mwSize        dim[3],
                          double        **ii,
                          double        **jj,
                          double        **nn,
                          double        **pp)
{
    int            cri = 0, nxti = 0;
    int            sl=0, r=0, c=0;
    double         cr = 0.0, nxt = 0.0;
    unsigned int   i = 0;
    size_t         size = 10000;
    
    (*ii) = mxCalloc(size,sizeof(double));
    (*jj) = mxCalloc(size,sizeof(double));
    (*nn) = mxCalloc(size,sizeof(double));
    (*pp) = mxCalloc(size,sizeof(double));
    
    for (sl=0; sl<dim[2]; sl++)
    {
        for (c=0; c<dim[1]; c++)
        {
            for (r=0; r<dim[0]; r++)
            {
                if ( (cr = rima[cri = index(r,c,sl,dim)]) )
                {
                    if (i > (size-4))
                    {
                        size += 10000;
                        (*ii) = mxRealloc((*ii),size*sizeof(double));
                        (*jj) = mxRealloc((*jj),size*sizeof(double));
                        (*nn) = mxRealloc((*nn),size*sizeof(double));
                        (*pp) = mxRealloc((*pp),size*sizeof(double));
                    }
                    if (sl)
                    {
                        if ((nxt = rima[nxti = index(r,c,sl-1,dim)]) && (cr != nxt))
                        {
                            (*ii)[i] = MIN(cr,nxt);
                            (*jj)[i] = MAX(cr,nxt);
                            (*nn)[i] = 1.0;
                            if (cr < nxt) (*pp)[i] = pm[cri] - pm[nxti];
                            else (*pp)[i] = pm[nxti] - pm[cri];
                            i++;
                        }
                    }
                    if (c)
                    {
                        if ((nxt = rima[nxti = index(r,c-1,sl,dim)]) && (cr != nxt))
                        {
                            (*ii)[i] = MIN(cr,nxt);
                            (*jj)[i] = MAX(cr,nxt);
                            (*nn)[i] = 1.0;
                            if (cr < nxt) (*pp)[i] = pm[cri] - pm[nxti];
                            else (*pp)[i] = pm[nxti] - pm[cri];
                            i++;
                        }
                    }
                    if (r)
                    {
                        if ((nxt = rima[nxti = index(r-1,c,sl,dim)]) && (cr != nxt))
                        {
                            (*ii)[i] = MIN(cr,nxt);
                            (*jj)[i] = MAX(cr,nxt);
                            (*nn)[i] = 1.0;
                            if (cr < nxt) (*pp)[i] = pm[cri] - pm[nxti];
                            else (*pp)[i] = pm[nxti] - pm[cri];
                            i++;
                        }
                    }
                }
            }
        }
    }

   return(i);
}


/* Gateway function with error check. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   mwSize         ndim, pm_ndim;
   int            n, i;
   const mwSize   *cdim = NULL, *pm_cdim = NULL;
   mwSize         dim[3];
   mwSize         nnz = 0;
   double         *rima = NULL;
   double         *pm = NULL;
   double         *ii = NULL, *jj = NULL;
   double         *nn = NULL, *pp = NULL;
   double         *tmp = NULL;


   if (nrhs == 0) mexErrMsgTxt("usage: [i,j,n,p]=pm_create_connectogram_dtj(rima,pm)");
   if (nrhs != 2) mexErrMsgTxt("pm_create_connectogram_dtj: 2 input arguments required");
   if (nlhs != 4) mexErrMsgTxt("pm_create_connectogram_dtj: 4 output argument required");

   /* Get connected components map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("pm_bwlabel_dtj: rima must be numeric, real, full and double");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if ((ndim < 2) | (ndim > 3))
   {
      mexErrMsgTxt("pm_bwlabel_dtj: rima must be 2 or 3-dimensional");
   }
   cdim = mxGetDimensions(prhs[0]);
   rima = mxGetPr(prhs[0]);

   /* Get phase-map. */

   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("pm_bwlabel_dtj: pm must be numeric, real, full and double");
   }
   pm_ndim = mxGetNumberOfDimensions(prhs[1]);
   if (pm_ndim != ndim)
   {
      mexErrMsgTxt("pm_bwlabel_dtj: rima and pm must have same dimensionality");
   }
   pm_cdim = mxGetDimensions(prhs[1]);
   for (i=0; i<ndim; i++)
   {
      if (cdim[i] != pm_cdim[i])
      {
         mexErrMsgTxt("pm_bwlabel_dtj: rima and pm must have same size");
      }
   }
   pm = mxGetPr(prhs[1]);

   /* Fix dimensions to allow for 2D and 3D data. */

   dim[0]=cdim[0]; dim[1]=cdim[1];
   if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 
   for (i=0, n=1; i<ndim; i++)
   {
      n *= dim[i];
   }

   /* 
      Create ii, jj, and nn and pp vectors for subsequent
      use by the matlab sparse function such that
      N = sparse(ii,jj,nn,nnz) generates a matrix where
      each non-zero entry signifies the no. of voxels
      along the border of the corresponding regions.
   */

   nnz = make_vectors(rima,pm,dim,&ii,&jj,&nn,&pp);
 
   /* Allocate memory for output. */

   plhs[0] = mxCreateDoubleMatrix(nnz,1,mxREAL);
   tmp = mxGetPr(plhs[0]);
   memcpy(tmp,ii,nnz*sizeof(double));
   plhs[1] = mxCreateDoubleMatrix(nnz,1,mxREAL);
   tmp = mxGetPr(plhs[1]);
   memcpy(tmp,jj,nnz*sizeof(double));
   plhs[2] = mxCreateDoubleMatrix(nnz,1,mxREAL);
   tmp = mxGetPr(plhs[2]);
   memcpy(tmp,nn,nnz*sizeof(double));
   plhs[3] = mxCreateDoubleMatrix(nnz,1,mxREAL);
   tmp = mxGetPr(plhs[3]);
   memcpy(tmp,pp,nnz*sizeof(double));

   /* Clean up a bit. */

   mxFree(ii); mxFree(jj); mxFree(nn); mxFree(pp);
   
}
