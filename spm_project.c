#ifndef lint
static char sccsid[]="%W% %E%";
#endif
 
/*

spm_project.c
% forms maximium intensity projections - a compiled routine
% FORMAT spm_project(X,L,SPM,V)
% X	-	a matrix of voxel values
% L	- 	a matrix of locations in Talairach et Tournoux (1988) space
% SPM	-	matrix for maximum intensity projection
% V     -       {1 x 6} vector of image and voxel sizes [DIM VOX] in mm
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

/* Input Arguments */
#define	V	prhs[0]
#define	L	prhs[1]
#define	SPM	prhs[2]
#define	DIM	prhs[3]

#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))


#ifdef __STDC__
void mexFunction(
	int		nlhs,
	Matrix	*plhs[],
	int		nrhs,
	Matrix	*prhs[]
	)
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
    double		*spm,*l,*v,*dim;
    unsigned int 	m,n,i,j,k;
    int			x,y,z,xdim,ydim,zdim;
    double		q;

    if (nrhs != 4 || nlhs > 0)
         mexErrMsgTxt("Inappropriate usage.");

    n    = mxGetN(V);
    m    = mxGetM(SPM);

    /* Assign pointers to the parameters */
    spm  = mxGetPr(SPM);
    l    = mxGetPr(L);
    v    = mxGetPr(V);
    dim  = mxGetPr(DIM);

    xdim = (int) (dim[3] + 0.5);
    ydim = (int) (dim[4] + 0.5);
    zdim = (int) (dim[5] + 0.5);

    /* go though point list */
    for (i = 0; i < n; i++) {
	x = (int) l[i*3 + 0];
	y = (int) l[i*3 + 1];
	z = (int) l[i*3 + 2];

	/* transverse */
	q = max(v[i], spm[(124 + y) + (104 - x)*m]);
	for (j = 0; j < ydim; j++) {
		for (k = 0; k < xdim; k++) {
		    	spm[124 + j + y + (104 + k - x)*m] = q;
		}
	}

	/* sagittal */
	q = max(v[i], spm[(124 + y) + (240 + z)*m]);
	for (j = 0; j < ydim; j++) {
		for (k = 0; k < zdim; k++) {
		    	spm[124 + j + y + (238 + k + z)*m] = q;
		}
	}

	/* coronal */
	q = max(v[i], spm[(276 + x) + (240 + z)*m]);
	for (j = 0; j < xdim; j++) {
		for (k = 0; k < zdim; k++) {
		    	spm[276 + j + x + (238 + k + z)*m] = q;
		}
	}

    }

}
