#include "mex.h"
#include <math.h>
#include <string.h>
#include <float.h>

#ifndef MAX
#define MAX(a,b)     (((a)>(b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)     (((a)<(b)) ? (a) : (b))
#endif

#ifndef INDX
#define INDX(x,y,z,dim) (((z)-1)*dim[0]*dim[1]+((y)-1)*dim[0]+((x)-1))
#endif

/* Functions doing the job. */

void invert_phasemap_dtj_idim1(mwSize *dim, double *pm, double *ipm)
{
    mwIndex x, y, z, i, oi;
    
    for (z=1; z<=dim[2]; z++)
    {
        for (y=1; y<=dim[1]; y++)
        {
            oi = 1;
            for (x=1; x<=dim[0]; x++)
            {
                for (i=oi; i<=dim[0] && (pm[INDX(i,y,z,dim)]+i)<x; i++) ; /* N.B. */
                if (i>1 && i<=dim[0])
                {
                    ipm[INDX(x,y,z,dim)] = ((double) i) - 1.0 - ((double) x) +
                            (((double) x)-pm[INDX(i-1,y,z,dim)]-((double) (i-1))) /
                            (pm[INDX(i,y,z,dim)]-pm[INDX(i-1,y,z,dim)]+1.0);
                }
                else
                {
                    ipm[INDX(x,y,z,dim)] = DBL_MAX;
                }
                oi = MAX(1,i-1);
            }
            
            for (i=1; i<=dim[0] && ipm[INDX(i,y,z,dim)]==DBL_MAX; i++) ; /* N.B. */
            for (x=i-1; x>0; x--)
            {
                ipm[INDX(x,y,z,dim)] = ipm[INDX(i,y,z,dim)];
            }
            
            for (i=dim[0]; i>0 && ipm[INDX(i,y,z,dim)]==DBL_MAX; i--) ; /* N.B. */
            for (x=i+1; x<=dim[0]; x++)
            {
                ipm[INDX(x,y,z,dim)] = ipm[INDX(i,y,z,dim)];
            }
        }
    }
}


void invert_phasemap_dtj_idim2(mwSize *dim, double *pm, double *ipm)
{
    mwIndex x, y, z, i, oi;
    
    for (z=1; z<=dim[2]; z++)
    {
        for (x=1; x<=dim[0]; x++)
        {
            oi = 1;
            for (y=1; y<=dim[1]; y++)
            {
                for (i=oi; i<=dim[1] && (pm[INDX(x,i,z,dim)]+i)<y; i++) ; /* N.B. */
                if (i>1 && i<=dim[1])
                {
                    ipm[INDX(x,y,z,dim)] = ((double) i) - 1.0 - ((double) y) +
                            (((double) y)-pm[INDX(x,i-1,z,dim)]-((double) (i-1))) /
                            (pm[INDX(x,i,z,dim)]-pm[INDX(x,i-1,z,dim)]+1.0);
                }
                else
                {
                    ipm[INDX(x,y,z,dim)] = DBL_MAX;
                }
                oi = MAX(1,i-1);
            }
            
            for (i=1; i<=dim[1] && ipm[INDX(x,i,z,dim)]==DBL_MAX; i++) ; /* N.B. */
            for (y=i-1; y>0; y--)
            {
                ipm[INDX(x,y,z,dim)] = ipm[INDX(x,i,z,dim)];
            }
            
            for (i=dim[1]; i>0 && ipm[INDX(x,i,z,dim)]==DBL_MAX; i--) ; /* N.B. */
            for (y=i+1; y<=dim[1]; y++)
            {
                ipm[INDX(x,y,z,dim)] = ipm[INDX(x,i,z,dim)];
            }
        }
    }
}


void invert_phasemap_dtj_idim3(mwSize *dim, double *pm, double *ipm)
{
    mwIndex x, y, z, i, oi;
    
    for (y=1; y<=dim[1]; y++)
    {
        for (x=1; x<=dim[0]; x++)
        {
            oi = 1;
            for (z=1; z<=dim[2]; z++)
            {
                for (i=oi; i<=dim[2] && (pm[INDX(x,y,i,dim)]+i)<z; i++) ; /* N.B. */
                if (i>1 && i<=dim[2])
                {
                    ipm[INDX(x,y,z,dim)] = ((double) i) - 1.0 - ((double) z) +
                            (((double) z)-pm[INDX(x,y,i-1,dim)]-((double) (i-1))) /
                            (pm[INDX(x,y,i,dim)]-pm[INDX(x,y,i-1,dim)]+1.0);
                }
                else
                {
                    ipm[INDX(x,y,z,dim)] = DBL_MAX;
                }
                oi = MAX(1,i-1);
            }
            
            for (i=1; i<=dim[2] && ipm[INDX(x,y,i,dim)]==DBL_MAX; i++) ; /* N.B. */
            for (z=i-1; z>0; z--)
            {
                ipm[INDX(x,y,z,dim)] = ipm[INDX(x,y,i,dim)];
            }
            
            for (i=dim[2]; i>0 && ipm[INDX(x,y,i,dim)]==DBL_MAX; i--) ; /* N.B. */
            for (z=i+1; z<=dim[2]; z++)
            {
                ipm[INDX(x,y,z,dim)] = ipm[INDX(x,y,i,dim)];
            }
        }
    }
}


/* Function doing the job. */

void invert_phasemap_dtj(mwSize *dim, double *pm, double *ipm)
{
    mwIndex i, oi;
    mwIndex f; /* Index in frequency encode direction. */
    mwIndex p; /* Index in phase encode direction. */
    mwIndex s; /* Index in slice selection direction. */
    
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwIndex i;
    mwSize dim[3] = {0, 0, 0};
    mwSize ndim, idim;
    
    if (nrhs == 0) mexErrMsgTxt("usage: ipm = pm_invert_phasemap_dtj(pm,idim)");
    if (nrhs != 2) mexErrMsgTxt("pm_invert_phasemap: 2 input argument required");
    if (nlhs != 1) mexErrMsgTxt("pm_invert_phasemap: 1 output arguments required");
    
    /* Get phasemap (deformation field). */
    
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("invert_phasemap_dtj: pm must be numeric, real, full and double");
    }
    
    if ((ndim = mxGetNumberOfDimensions(prhs[0])) > 3)
    {
        mexErrMsgTxt("invert_phasemap_dtj: pm must be 2 or 3 dimensional");
    }
    for (i=0;i<ndim;i++)
    {
        dim[i] = mxGetDimensions(prhs[0])[i];
    }
    for (i=ndim; i<3; i++)
    {
        dim[i] = 1;
    }
    
    /* Get dimension along which to perform inversion. */
    
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("invert_phasemap_dtj: idim must be numeric, real, full and double");
    }
    if (mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("invert_phasemap_dtj: idim must be scalar");
    }
    idim = (mwSize) mxGetScalar(prhs[1]);
    if (idim < 1 || idim > 3)
    {
        mexErrMsgTxt("invert_phasemap_dtj: idim must be in range 1 - 3");
    }
    
    /* Allocate memory for output. */
    
    plhs[0] = mxCreateNumericArray(ndim,dim,mxDOUBLE_CLASS,mxREAL);
    
    if (idim == 1)
    {
        invert_phasemap_dtj_idim1(dim,mxGetPr(prhs[0]),mxGetPr(plhs[0]));
    }
    else if (idim == 2)
    {
        invert_phasemap_dtj_idim2(dim,mxGetPr(prhs[0]),mxGetPr(plhs[0]));
    }
    else if (idim == 3)
    {
        invert_phasemap_dtj_idim3(dim,mxGetPr(prhs[0]),mxGetPr(plhs[0]));
    }
    
}
