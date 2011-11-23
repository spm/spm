/*
 * $Id: spm_dilate_erode.c 4569 2011-11-23 16:11:30Z guillaume $
 * Jesper Andersson
 */

#include <string.h>
#include "mex.h"

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

/* Utility function that returns index into */
/* 1D array with range checking.            */
 
mwSignedIndex get_index(mwSignedIndex i,
                        mwSignedIndex j,
                        mwSignedIndex k,
                        mwSize        dim[3])
{
    if ((i<0) || (i>(dim[0]-1)) || (j<0) || (j>(dim[1]-1)) || (k<0) || (k>(dim[2]-1)))
        return(-1);
    else
        return(k*dim[0]*dim[1]+j*dim[0]+i);
}

/* Function doing the job. */

void do_it(double  *iima,
           mwSize  dim[3],
           double  *krnl,
           mwSize  kdim[3],
           int     dilate,
           double  *oima)
{
    mwIndex  i=0, j=0, k=0;
    mwIndex  ki=0, kj=0, kk=0;
    mwIndex  ndx=0;
    mwSignedIndex kndx=0;
    double   kv=0.0;
  
    for (i=0; i<dim[0]; i++)
    {
        for (j=0; j<dim[1]; j++)
        {
            for (k=0; k<dim[2]; k++)
            {
                ndx = get_index(i,j,k,dim);
                for (ki=0; ki<kdim[0]; ki++)
                {
                    for (kj=0; kj<kdim[1]; kj++)
                    {
                        for (kk=0; kk<kdim[2]; kk++)
                        {
                            kndx = get_index(i-(kdim[0]/2)+ki,j-(kdim[1]/2)+kj,k-(kdim[2]/2)+kk,dim);
                            if (kndx > -1)
                            {
                                if (dilate)
                                {
                                    if ((kv=krnl[get_index(ki,kj,kk,kdim)]))
                                    {
                                        oima[ndx] = MAX(oima[ndx],kv*iima[kndx]);
                                    }
                                }
                                else /* erode */
                                {
                                    if ((kv=krnl[get_index(ki,kj,kk,kdim)]))
                                    {
                                        oima[ndx] = MIN(oima[ndx],kv*iima[kndx]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
   }
}


/* Gateway function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char          *fnc_str = NULL;
    mwSize        ndim=0, krn_ndim=0;
    int           n, i;
    int           dilate = 0;
    mwSize        buflen = 0;
    const mwSize  *cdim = NULL, *krn_cdim = NULL;
    mwSize        dim[3], kdim[3];
    double        *iima = NULL;
    double        *oima = NULL;
    double        *krnl = NULL;

    if (nrhs < 3) mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 3) mexErrMsgTxt("Too many input arguments.");
    if (nlhs < 1) mexErrMsgTxt("Not enough output arguments");
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");

    /* Get image */

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("spm_dilate_erode: ima must be numeric, real, full and double");
    }
    ndim = mxGetNumberOfDimensions(prhs[0]);
    if ((ndim < 2) || (ndim > 3))
    {
        mexErrMsgTxt("spm_dilate_erode: ima must be 2 or 3-dimensional");
    }
    cdim = mxGetDimensions(prhs[0]);
    iima = mxGetPr(prhs[0]);

    /* Fix dimensions to allow for 2D and 3D data */

    dim[0] = cdim[0]; dim[1] = cdim[1];
    if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 
    for (i=0, n=1; i<ndim; i++)
    {
        n *= dim[i];
    }

    /* Get kernel */

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("spm_dilate_erode: kernel must be numeric, real, full and double");
    }
    krn_ndim = mxGetNumberOfDimensions(prhs[1]);
    if (krn_ndim != mxGetNumberOfDimensions(prhs[0]))
    {
        mexErrMsgTxt("spm_dilate_erode: ima and kernel must have same dimensionality");
    }
    krn_cdim = mxGetDimensions(prhs[1]);
    krnl = mxGetPr(prhs[1]);
    kdim[0]=krn_cdim[0]; kdim[1]=krn_cdim[1];
    if (krn_ndim==2) {kdim[2]=1; krn_ndim=3;} else {kdim[2]=krn_cdim[2];} 

    /* Check if dilate or erode */

    if (!mxIsChar(prhs[2]) || (mxGetM(prhs[2]) != 1))
    {
        mexErrMsgTxt("spm_dilate_erode: function should be a string");
    }
   
    buflen = mxGetN(prhs[2])*sizeof(mxChar)+1;
    fnc_str = (char *) mxMalloc(buflen);
    mxGetString(prhs[2],fnc_str,buflen);

    if (strcmp(fnc_str,"dilate") && strcmp(fnc_str,"erode"))
    {
        mexErrMsgTxt("spm_dilate_erode: function must have value 'dilate' or 'erode'");
    }
    if (!strcmp(fnc_str,"dilate")) {dilate = 1;}
    else {dilate = 0;}
   
    /* Allocate and initialise output image */

    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                   mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
    oima = mxGetPr(plhs[0]);
    memcpy(oima,iima,n*sizeof(double));

    /*   if (dilate) {for (i=0; i<n; i++) oima[i]=-DBL_MAX;}
    else {for (i=0; i<n; i++) oima[i]=DBL_MAX;} */

    do_it(iima,dim,krnl,kdim,dilate,oima);

    mxFree(fnc_str);
}
