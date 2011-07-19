/*
 * $Id: spm_project.c 4399 2011-07-19 16:28:36Z guillaume $
 */
 
/*

spm_project.c
% forms maximium intensity projections - a compiled routine
% FORMAT spm_project(X,L,dims)
% X    - a matrix of voxel values
% L    - a matrix of locations in Talairach et Tournoux (1988) space
% dims - assorted dimensions
%        dims(1:3) - the sizes of the projected rectangles.
%        dims(4:5) - the dimensions of the mip image.
%____________________________________________________________________________
%
% spm_project 'fills in' a matrix (SPM) in the workspace to create
% a maximum intensity projection according to a point list of voxel
% values (V) and their locations (L) in the standard space described
% in the atlas of Talairach & Tournoux (1988).
%
% see also spm_mip.m

*/

#include <math.h>
#include <stdio.h>
#include "mex.h"
#define RINT(A) floor((A)+0.5)
#define MAX(A, B)   ((A) > (B) ? (A) : (B))
#define MIN(A, B)   ((A) < (B) ? (A) : (B))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *spm,*l,*v,*dim;
    int    m,m1,n,i,j,k, o;
    int    x,y,z,xdim,ydim,zdim;
    double q;
    double *DXYZ, *CXYZ;
    int    DX, DY, DZ, CX, CY, CZ;

    if ( !(nrhs == 3 || nrhs == 5) || nlhs > 1) mexErrMsgTxt("Incorrect usage.");

    for(k=0; k<nrhs; k++)
        if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
            mxIsSparse(prhs[k]) || !mxIsDouble(prhs[k]))
            mexErrMsgTxt("Arguments must be numeric, real, full and double.");

    /* The values */
    n    = mxGetNumberOfElements(prhs[0]);
    v    = mxGetPr(prhs[0]);

    /* The co-ordinates */
    if ((mxGetN(prhs[1]) != n) || (mxGetM(prhs[1]) != 3))
        mexErrMsgTxt("Incompatible size for locations matrix.");
    l    = mxGetPr(prhs[1]);

    /* Dimension information */
    if (mxGetNumberOfElements(prhs[2]) != 5)
        mexErrMsgTxt("Incompatible size for dimensions vector.");
    dim  = mxGetPr(prhs[2]);
    xdim = (int) (fabs(dim[0]) + 0.99);
    ydim = (int) (fabs(dim[1]) + 0.99);
    zdim = (int) (fabs(dim[2]) + 0.99);
    m    = (int) (dim[3]);
    m1   = (int) (dim[4]);

    if (nrhs == 3)
    {
      DX = 182;
      DY = 218;
      DZ = 182;
      CX =  89;
      CY = 125;
      CZ =  71;
    }
    else
    {
      if (mxGetNumberOfElements(prhs[3]) != 3)
        mexErrMsgTxt("Incompatible size for DXYZ vector.");
      DXYZ = mxGetPr(prhs[3]);
      DX = (int) (DXYZ[0]);
      DY = (int) (DXYZ[1]);
      DZ = (int) (DXYZ[2]);
      if (mxGetNumberOfElements(prhs[4]) != 3)
        mexErrMsgTxt("Incompatible size for CXYZ vector.");
      CXYZ = mxGetPr(prhs[4]);
      CX = (int) (CXYZ[0]);
      CY = (int) (CXYZ[1]);
      CZ = (int) (CXYZ[2]);
    }

    plhs[0] = mxCreateDoubleMatrix(m,m1,mxREAL);
    spm     = mxGetPr(plhs[0]);

    if (m == DY+DX && m1 == DZ+DX) /* MNI Space */
    {
        /* go though point list */
        for (i = 0; i < n; i++)
        {
            x = (int)RINT(l[i*3 + 0]) + CX;
            y = (int)RINT(l[i*3 + 1]) + CY;
            z = (int)RINT(l[i*3 + 2]) + CZ;

            if (2*CX-x-xdim/2>=0 && 2*CX-x+xdim/2<DX && y-ydim/2>=0 && y+ydim/2<DY) /* transverse */
            {
                q = v[i];
                for (j = -ydim/2; j <= ydim/2; j++)
                    for (k = -xdim/2; k <= xdim/2; k++)
                    {
                        o = j + y + (k + 2*CX-x)*m;
                        if (spm[o]<q) spm[o] = q;
                    }
            }

            if (z-zdim/2>=0 && z+zdim/2<DZ && y-ydim/2>=0 && y+ydim/2<DY) /* sagittal */
            {
                q = v[i];
                for (j = -ydim/2; j <= ydim/2; j++)
                    for (k = -zdim/2; k <= zdim/2; k++)
                    {
                        o = j + y + (DX + k + z)*m;
                        if (spm[o]<q) spm[o] = q;
                    }
            }

            if (x-xdim/2>=0 && x+xdim/2<DX && z-zdim/2>=0 && z+zdim/2<DZ) /* coronal */
            {
                q = v[i];
                for (j = -xdim/2; j <= xdim/2; j++)
                    for (k = -zdim/2; k <= zdim/2; k++)
                    {
                        o = DY + j + x + (DX + k + z)*m;
                        if (spm[o]<q) spm[o] = q;
                    }
            }
        }
    }
    else if (m == 360 && m1 == 352) /* old code for the old MIP matrix */
    {
    for (i = 0; i < n; i++) {
        x = (int) l[i*3 + 0];
        y = (int) l[i*3 + 1];
        z = (int) l[i*3 + 2];
    
        /* transverse */
        q = MAX(v[i], spm[(124 + y) + (104 - x)*m]);
        for (j = 0; j < ydim; j++) {
            for (k = 0; k < xdim; k++) {
                spm[124 + j + y + (104 + k - x)*m] = q;
            }
        }
    
        /* sagittal */
        q = MAX(v[i], spm[(124 + y) + (240 + z)*m]);
        for (j = 0; j < ydim; j++) {
            for (k = 0; k < zdim; k++) {
                spm[124 + j + y + (238 + k + z)*m] = q;
            }
        }
    
        /* coronal */
        q = MAX(v[i], spm[(276 + x) + (240 + z)*m]);
        for (j = 0; j < xdim; j++) {
            for (k = 0; k < zdim; k++) {
                spm[276 + j + x + (238 + k + z)*m] = q;
            }
        }
    }
    }
    else 
        mexErrMsgTxt("Wrong sized MIP matrix");
}
