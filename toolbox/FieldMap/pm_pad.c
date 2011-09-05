#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <limits.h>

void pad(double         *pm,
         double         *wmap,
         mwSize         dim[3],
         double         *krnl,
         mwSize         kdim[3],
         double         *opm,
         double         *owmap);

int getindex(mwSignedIndex   i,
             mwSignedIndex   j,
             mwSignedIndex   k,
             mwSize          dim[3]);


void pad(double         *pm,
         double         *wmap,
         mwSize         dim[3],
         double         *krnl,
         mwSize         kdim[3],
         double         *opm,
         double         *owmap)
{
    mwIndex        i=0, j=0, k=0;
    mwIndex        ki=0, kj=0, kk=0;
    int            n = 0;
    int            ndx=0, kndx=0;
    double         p = 0.0;
    
    for (i=0; i<dim[0]; i++)
    {
        for (j=0; j<dim[1]; j++)
        {
            for (k=0; k<dim[2]; k++)
            {
                if (!wmap[ndx=getindex(i,j,k,dim)])
                {
                    n = 0;
                    p = 0.0;
                    for (ki=0; ki<kdim[0]; ki++)
                    {
                        for (kj=0; kj<kdim[1]; kj++)
                        {
                            for (kk=0; kk<kdim[2]; kk++)
                            {
                                kndx = getindex(i-(kdim[0]/2)+ki,j-(kdim[1]/2)+kj,k-(kdim[2]/2)+kk,dim);
                                if (kndx > -1)
                                {
                                    if (wmap[kndx])
                                    {
                                        p += krnl[getindex(ki,kj,kk,kdim)] * pm[kndx];
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                    if (n)
                    {
                        opm[ndx] = p/n;
                        owmap[ndx] = 1;
                    }
                    else
                    {
                        opm[ndx]=pm[ndx];
                        owmap[ndx]=wmap[ndx];
                    }
                }
                else
                {
                    opm[ndx]=pm[ndx];
                    owmap[ndx]=wmap[ndx];
                }
            }
        }
    }
    return;
}

/* Utility function that returns index into */
/* 1D array with range checking.            */
 
int getindex(mwSignedIndex  i,
             mwSignedIndex  j,
             mwSignedIndex  k,
             mwSize         dim[3])
{
   if ((i<0) | (i>(dim[0]-1)) | 
       (j<0) | (j>(dim[1]-1)) | 
       (k<0) | (k>(dim[2]-1)))
       return(-1);
   else
       return(k*dim[0]*dim[1]+j*dim[0]+i);
}


/* Gateway function with error check. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   mwSize         ndim, wmap_ndim, krn_ndim;
   mwIndex        n, i;
   const mwSize   *cdim = NULL, *wmap_cdim = NULL, *krn_cdim = NULL;
   mwSize         dim[3], kdim[3];
   double         *pm = NULL;
   double         *wmap = NULL;
   double         *opm = NULL;
   double         *owmap = NULL;
   double         *krnl = NULL;


   if (nrhs == 0) mexErrMsgTxt("usage: [pm,wmap] = pm_pad(pm,wmap,kernel)");
   if (nrhs != 3) mexErrMsgTxt("pm_pad: 3 input arguments required");
   if (nlhs != 2) mexErrMsgTxt("pm_pad: 2 output argument required");

   /* Get phase map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("pm_pad: pm must be numeric, real, full and double");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if ((ndim < 2) | (ndim > 3))
   {
      mexErrMsgTxt("pm_pad: pm must be 2 or 3-dimensional");
   }
   cdim = mxGetDimensions(prhs[0]);
   pm = mxGetPr(prhs[0]);

   /* Get wrap-map. */

   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("pm_pad: wmap must be numeric, real, full and double");
   }
   wmap_ndim = mxGetNumberOfDimensions(prhs[1]);
   if (wmap_ndim != ndim)
   {
      mexErrMsgTxt("pm_pad: pm and wmap must have same dimensionality");
   }
   wmap_cdim = mxGetDimensions(prhs[1]);
   for (i=0; i<ndim; i++)
   {
      if (cdim[i] != wmap_cdim[i])
      {
         mexErrMsgTxt("pm_pad: pm and wmap must have same size");
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
      mexErrMsgTxt("pm_pad: kernel must be numeric, real, full and double");
   }
   krn_ndim = mxGetNumberOfDimensions(prhs[2]);
   /*if (krn_ndim != ndim)
   {
      mexErrMsgTxt("pm_pad: pm and kernel must have same dimensionality");
      }*/
   krn_cdim = mxGetDimensions(prhs[2]);
   krnl = mxGetPr(prhs[2]);
   kdim[0]=krn_cdim[0]; kdim[1]=krn_cdim[1];
   if (krn_ndim==2) {kdim[2]=1; krn_ndim=3;} else {kdim[2]=krn_cdim[2];} 
   if (krn_ndim != ndim)
   {
      mexErrMsgTxt("pm_pad: pm and kernel must have same dimensionality");
      }
   /* Allocate padded output phasemap end new wmap. */

   plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                  mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
   opm = mxGetPr(plhs[0]);

   plhs[1] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                  mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
   owmap = mxGetPr(plhs[1]);

   pad(pm,wmap,dim,krnl,kdim,opm,owmap);
   
}
