#ifndef lint
static char sccsid[]="@(#)spm_max.c	2.4 JB Poline 99/05/18";
#endif

/*  spm_max.c JB Poline 10/11/94 

% Return sizes, maxima and locations of local excursion {X > u} sets 
% FORMAT [N Z M A] = spm_lmax(X,L)
%
% X     - values of 3-D field
% L     - locations [x y x]' {in vox}
%
% N     - size of region {in voxels)
% Z     - Z values of maxima
% M     - location of maxima {in vox}
% A     - region number
%____________________________________________________________________________
%
% spm_max characterizes a point list of voxel values (X) and their locations
% (L) in terms of edge, face and vertex connected subsets, returning a maxima-
% orientated list:  The value of the ith maximum is Z(i) and its location
% is given by M(:,i). A(i) identifies the ith maximum with a region. Region
% A(i) contains N(i) voxels.
% 
*/

#include <math.h>
#include "mex.h"

/* Input Output Arguments */
#define	SPM	prhs[0]
#define	XYZ	prhs[1]

#define	NUM	plhs[0]
#define	Z_MAX	plhs[1]
#define	LOC	plhs[2]
#define	IND	plhs[3]


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

void find_loc_max_strict_c18(float *vol, int dx, int dy, int dz,
	int minx, int miny, int minz, double *max, double *loc, int *nb_lm)
{
	int	i,j,k, dxdy = dx*dy;
	float	*Uv, *Lv, *Cv, *Uvp,
		*Lvp, *Cvp, *Uvl, *Lvl, *Cvl;


	/*---------------------- */
     	/* find the local maxima */
	
	*nb_lm = 0; 			/* nb_lm : nb of loc max */
	for( k=1; k<dz-1; k++)
	{
	 Cvp = vol + k*dxdy;
	 Lvp = Cvp - dxdy;
	 Uvp = Cvp + dxdy;
	    for( j=1; j<dy-1; j++)
	    {
	     Cvl = Cvp + j*dx;
	     Lvl = Lvp + j*dx;
	     Uvl = Uvp + j*dx;
	     
		for(i=1; i<dx-1; i++)
		{
		 Cv = Cvl + i;
	     	 Lv = Lvl + i;
	     	 Uv = Uvl + i;

		/* 18 connectivity */

		 if( 	( *(Cv + 1) < *Cv ) &&   /* East */
			( *(Cv - 1) < *Cv ) &&   /* West  */
			( *(Cv + dx) < *Cv ) &&  /* south */
			( *(Cv - dx) < *Cv ) &&  /* north */

			( *(Cv + dx + 1) < *Cv ) &&  /* SE */
			( *(Cv + dx - 1) < *Cv ) &&  /* SW */
			( *(Cv - dx + 1) < *Cv ) &&  /* NE */
			( *(Cv - dx - 1) < *Cv ) &&  /* NW */


			( *Lv < *Cv ) &&	 /* lower */ 
			( *(Lv + 1) < *Cv ) &&	 /* LE */ 
			( *(Lv - 1) < *Cv ) &&	 /* LW */ 
			( *(Lv + dx) < *Cv ) &&	 /* LS */ 
			( *(Lv - dx) < *Cv ) &&	 /* LN */ 

 			( *Uv < *Cv ) && 	 /* upper */
			( *(Uv + 1) < *Cv ) &&	 /* UE */ 
			( *(Uv - 1) < *Cv ) &&	 /* UW */ 
			( *(Uv + dx) < *Cv ) &&	 /* US */ 
			( *(Uv - dx) < *Cv ) 	 /* UN */ 
		   )

		  {
		   loc[(3* (*nb_lm) + 0)] = (i + minx -1);
	    	   loc[(3* (*nb_lm) + 1)] = (j + miny -1);
	    	   loc[(3* (*nb_lm) + 2)] = (k + minz -1);
	    	   max[ *nb_lm ]         = (double) *Cv;
	    	   (*nb_lm)++;		 		


		  } /* if */
		}
	    }
	}

} /* find_loc_max_strict_c18 */



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double		*X,*L,*loc,*max;
    float		*vol, fl_tmp;
    int 		n,i,j,d;
    int			max_x, max_y, max_z, min_x, min_y, min_z,
			dx, dy, dz, dxdy, *num, *ind, index, nb_lm;
    int			*x,*y,*z;

    if ((nrhs != 2 ) || (nlhs > 4))
		mexErrMsgTxt("Inappropriate usage.");

    /* Assign pointers to the parameters */
    X    = mxGetPr(SPM);
    L    = mxGetPr(XYZ);


    i    = mxGetM(XYZ);
    if(i != 3) mexErrMsgTxt("XYZ M-dimension should be 3 ");

    n    = mxGetN(XYZ);
    i    = mxGetM(SPM);
    j    = mxGetN(SPM);

	
    if(!n) mexErrMsgTxt("no points in the volume in spm_max");
    if((i*j) != n) mexErrMsgTxt(" XYZ SPM dimension non compatible ");


    /* set up point list in x y z */
    x    = (int *) mxCalloc (n,sizeof( int ));
    y    = (int *) mxCalloc (n,sizeof( int ));
    z    = (int *) mxCalloc (n,sizeof( int ));


    max_x =  (int) floor( L[0] + 0.5);
    max_y =  (int) floor( L[1] + 0.5);
    max_z =  (int) floor( L[2] + 0.5);
    min_x =  (int) floor( L[0] + 0.5);
    min_y =  (int) floor( L[1] + 0.5);
    min_z =  (int) floor( L[2] + 0.5);

    for (i = 0, d = 0; i < 3*n; i = i + 3, d++) {
	x[d] = (int) floor(L[i + 0] + 0.5);
	if (x[d] > max_x) max_x = x[d];
	if (x[d] < min_x) min_x = x[d];
	y[d] = (int) floor(L[i + 1] + 0.5);
	if (y[d] > max_y) max_y = y[d];
	if (y[d] < min_y) min_y = y[d];
	z[d] = (int) floor(L[i + 2] + 0.5);
	if (z[d] > max_z) max_z = z[d];
	if (z[d] < min_z) min_z = z[d];	
    }


    /* create a volume of the points plus a bounding box */
	dx = max_x - min_x + 3;
	dy = max_y - min_y + 3;
	dz = max_z - min_z + 3;

	dxdy = dx*dy;

	vol  = (float *) mxCalloc((dx*dy*dz), sizeof(*vol));
	if(!vol) mexErrMsgTxt("\n memory alloc pb in loc_max ");

	/* put the points in the volume */
	for( i=0; i<d; i++) {
	  *(vol + (z[i] - min_z +1)*dxdy
		+ (y[i] - min_y +1)*dx 
		+  x[i] - min_x +1) = (float) X[i];
	}
	

    	max  = (double *) mxCalloc (n,sizeof( double ));
    	loc  = (double *) mxCalloc (n*3,sizeof( double ));
    	num  = (int *)    mxCalloc (n,sizeof( int ));
   	ind  = (int *)    mxCalloc (n,sizeof( int ));

	if (!num || !loc || !max || !ind)
		mexErrMsgTxt("\n memory alloc pb in loc_max ");


	find_loc_max_strict_c18(vol,dx,dy,dz, min_x,min_y,min_z, max,loc, &nb_lm);


/*--------------------------------*/
/* find the size of the local max */
	
    /* global settings */
    	pos3.t      = EPS2;	/* threshold + epsilon */
    	pos3.dx     = dx;
    	pos3.dy     = dy;
    	pos3.dz     = dz;

    index = 1;
    for (i = 0; i < nb_lm; i++) {

	pos3.x = (int) floor( loc[i*3 + 0] - min_x + 1+ 0.5);
	pos3.y = (int) floor( loc[i*3 + 1] - min_y + 1+ 0.5);
	pos3.z = (int) floor( loc[i*3 + 2] - min_z + 1+ 0.5);
	pos3.size  = 0;
	pos3.label = (float) -(i+1) ;
	
	/*--------------------------------*/
	/* trick : replace the voxels of the volume which	*/
	/* have been passed through by -(i+1);			*/

	if( (fl_tmp = *(vol + pos3.z*dxdy + pos3.y*dx + pos3.x)) < 0 ) {

	  /* point already belong to a cluster */
	  /* put the size of this cluster      */
	  j =  (int) floor( fl_tmp+ 0.5);
	  num[i] = num[-j-1];
	  ind[i] = ind[-j-1];
	}
	else {
	  /* compute cluster size */
	  pos3.ptr_flt = (vol + pos3.z*dxdy + pos3.y*dx + pos3.x);
	  comp_c_c18();
	  num[i] = pos3.size;
	  ind[i] = index;
	  index++;
	}
	

    } /* for */

    mxFree(vol);

    Z_MAX   = mxCreateDoubleMatrix(1, nb_lm, mxREAL);
    LOC   = mxCreateDoubleMatrix(3, nb_lm, mxREAL);
    NUM   = mxCreateDoubleMatrix(1, nb_lm, mxREAL);
    IND   = mxCreateDoubleMatrix(1, nb_lm, mxREAL);


    for (i = 0; i < nb_lm; i++) {
	    mxGetPr(NUM)[i] = (double) num[i];
    }
    for (i = 0; i < nb_lm; i++) {
	    mxGetPr(Z_MAX)[i] = (double) max[i];
    }
    for (i = 0; i < 3*nb_lm; i++) {
	    mxGetPr(LOC)[i] = (double) loc[i];
    }
    for (i = 0; i < nb_lm; i++) {
	    mxGetPr(IND)[i] = (double) ind[i];
    }
    mxFree(num);
    mxFree(loc);
    mxFree(max);
    mxFree(ind);
}
