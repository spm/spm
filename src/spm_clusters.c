#ifndef lint
static char sccsid[]="@(#)spm_clusters.c	2.2 JB Poline 99/03/03";
#endif

/*  spm_plateau.c JB Poline 10/11/94 

% Return cluster index for a point list
% FORMAT [A] = spm_cluster2(L)
% L     - locations [x y x]' {in vox}
%
% A     - cluster index or region number
%____________________________________________________________________________
%
% spm_clusters characterizes a point list of voxel values (X) and their 
% locations (L) in terms of edge, face and vertex connected subsets, returning a
% list of indices in A, such that X(i) belongs to cluster A(i) {using an 18 
% connectivity scheme)
% 
*/

#include <math.h>
#include "mex.h"

/* Input Output Arguments */
#define	XYZ	prhs[0]

#define	IND	plhs[0]


#include "connex.h"	

	/* where 
	
	typedef struct cc_position_3d_str
		{
		float	t, sum, label,  *ptr_flt;
		int	x, y, z, dx, dy, dz, size;
		}
		cc_position_3d;
	
	
	static	cc_position_3d	pos3;
	
	is defined */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double		*L;
    float		*vol, fl_tmp;
    int 		i,j,d,r;
    int			max_x, max_y, max_z, min_x, min_y, min_z,
			 dx, dy, dz, dxdy, *ind, index;
    int			*x,*y,*z;


    if ((nrhs != 1) || (nlhs > 1))
		mexErrMsgTxt("Inappropriate usage.");

    /* Assign pointers to the parameters */
    L    = mxGetPr(XYZ);
    i    = mxGetM(XYZ);
    if(i != 3) mexErrMsgTxt("  XYZ M-dimension should be 3 ");
    r    = mxGetN(XYZ);
    if(!r) mexErrMsgTxt("no points in the volume in spm_max");

    /* set up point list in x y z */
    x    = (int *) mxCalloc (r,sizeof( int ));
    y    = (int *) mxCalloc (r,sizeof( int ));
    z    = (int *) mxCalloc (r,sizeof( int ));

    ind  = (int *)    mxCalloc (r,sizeof( int ));

    if(!x || !y || !z || !ind) mexErrMsgTxt("\n memory alloc pb in loc_max ");

    max_x =  (int) floor(L[0] + 0.5);
    max_y =  (int) floor(L[1] + 0.5);
    max_z =  (int) floor(L[2] + 0.5);
    min_x =  (int) floor(L[0] + 0.5);
    min_y =  (int) floor(L[1] + 0.5);
    min_z =  (int) floor(L[2] + 0.5);

    d	 = 0;
    for (i = 0; i < 3*r; i = i + 3) {
	x[d] = (int) floor(L[i + 0] + 0.5);
	if (x[d] > max_x) max_x = x[d];
	if (x[d] < min_x) min_x = x[d];
	y[d] = (int) floor(L[i + 1] + 0.5);
	if (y[d] > max_y) max_y = y[d];
	if (y[d] < min_y) min_y = y[d];
	z[d] = (int) floor(L[i + 2] + 0.5);
	if (z[d] > max_z) max_z = z[d];
	if (z[d] < min_z) min_z = z[d];
	d++;
    }


    /* create a volume of the points plus a bounding box */
	dx = max_x - min_x + 3;
	dy = max_y - min_y + 3;
	dz = max_z - min_z + 3;

	dxdy = dx*dy;

	vol  = (float *) mxCalloc((dx*dy*dz), sizeof(*vol));
	if(!vol) mexErrMsgTxt("\n memory alloc pb in loc_max ");

	/* put the points in the volume */
	for( i=0; i<r; i++) {
	  *(vol + (z[i] - min_z + 1)*dxdy
		+ (y[i] - min_y + 1)*dx 
		+  x[i] - min_x + 1) = (float) 1.0;
	}
	
	
    /* global settings */
    	pos3.t      = EPS2;	/* threshold + epsilon */
    	pos3.dx     = dx;
    	pos3.dy     = dy;
    	pos3.dz     = dz;

    index = 1;
    for (i = 0; i < r; i++) {

	pos3.x = (int) floor( L[i*3 + 0] - min_x + 1 + 0.5);
	pos3.y = (int) floor( L[i*3 + 1] - min_y + 1 + 0.5);
	pos3.z = (int) floor( L[i*3 + 2] - min_z + 1 + 0.5);
	pos3.size  = 0;
	pos3.label = (float) -(i+1) ;
	
	/*--------------------------------*/
	/* trick : replace the voxels of the volume which	*/
	/* have been passed through by -(i+1);			*/

	if( (fl_tmp = *(vol + pos3.z*dxdy + pos3.y*dx + pos3.x)) < EPS2 ) {

	  /* point already belong to a cluster */
	  j =  (int) floor(fl_tmp+0.5);
	  ind[i] = ind[-j-1];
	}

	else {

	  /* compute cluster size */
	  pos3.ptr_flt = (vol + pos3.z*dxdy + pos3.y*dx + pos3.x);
	  comp_c_c18();

	  ind[i] = index;
	  index++;
	}
	

    } /* for */

    mxFree(vol);

    IND   = mxCreateDoubleMatrix(1, r, mxREAL);

    for (i = 0; i < r; i++) {
	    mxGetPr(IND)[i] = (double) ind[i];
    }
    mxFree(ind);


}
