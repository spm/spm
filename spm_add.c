#ifndef lint
static char sccsid[]="%W% Jean-Baptiste Poline, John Ashburner %E%";
#endif
 
/*
spm_mean.c
% averages a series of images - a compiled routine
% FORMAT s = spm_mean(n,TYPE,Q,P,S)
% n	-	number of voxels per image
% TYPE	- 	data type (see spm_type.m & volume.h)
% Q	-	filename for averaged image
% P	-	matrix of filenames to be averaged (rowwise strings)
% S	-	optional vector of scalefactors / weights
% s	-	the scalefactor which should be assigned to Q
%_______________________________________________________________________
%
% spm_mean computes a weighted sum of a set of image files to produce
% an mean image that is written to a named file (Q). No headers are
% read or written, in particular scalefactors are ignored. The image is
% written in the same type as the input images.
%
% The weights (S) default to 1/size(P,1) - resulting an average of the
% image files being written to Q. For a "proper" average including
% scalefactors, the weights s should be specified as sf/length(sf) for
% sf a vector of scalefactors for the image files specified in P.
%
% See also: spm_average - for details on how to usefully use this function
%
%_______________________________________________________________________

*/

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "mex.h"
#include "volume.h"


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
	MAPTYPE **maps;
	double *sptr, *scales, *image, scale;
	short **dptr;
	int ni, nj, nk, i, j, k, fo;
	static double mat[] = {1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1};
	char label[1024];

	if (nrhs != 2 || nlhs > 1)
		mexErrMsgTxt("Inappropriate usage.");

	maps = get_maps(prhs[0], &ni);


	for(i=1; i<ni; i++)
	{
		if (	maps[i]->xdim != maps[0]->xdim ||
			maps[i]->ydim != maps[0]->ydim ||
			maps[i]->zdim != maps[0]->zdim)
			{
				mexErrMsgTxt("Incompatible image dimensions.");
			}

		if (	maps[i]->xpixdim != maps[0]->xpixdim ||
			maps[i]->ypixdim != maps[0]->ypixdim ||
			maps[i]->zpixdim != maps[0]->zpixdim)
			{
				mexErrMsgTxt("Incompatible voxel sizes.");
			}
	}

	mxGetString(prhs[1],label,mxGetN(prhs[1])+1);
	fo  = open(label, O_RDWR|O_CREAT, 0644);
	if (fo == -1)
		mexErrMsgTxt("Cant open output file.");

	nj = maps[0]->zdim;
	nk = maps[0]->xdim*maps[0]->ydim;

	scales = (double *)mxCalloc(nj, sizeof(double));
	sptr   = (double *)mxCalloc(nk, sizeof(double));
	image  = (double *)mxCalloc(nk, sizeof(double));
	dptr   = (short **)mxCalloc(nj, sizeof(double));

	for(j=0; j<maps[0]->zdim; j++)
	{
		double mx, mn;
		mat[14] = j+1.0;

		for(k=0; k<nk; k++)
			sptr[k] = 0.0;

		for(i=0; i<ni; i++)
		{
			slice(mat, image, (int)maps[0]->xdim,(int)maps[0]->ydim, maps[i]->data,
				(int)maps[i]->xdim,(int)maps[i]->ydim,(int)maps[i]->zdim, 0, 0.0,
				maps[i]->scalefactor,0.0, maps[i]->datatype);
			for(k=0; k<nk; k++)
				sptr[k] += image[k];
		}

		/* Determine maximum and minimum */
		mx = -9e99;
		mn = 9e99;
		for(k=0; k<nk; k++)
		{
			if (sptr[k]>mx) mx=sptr[k];
			if (sptr[k]<mn) mn=sptr[k];
		}

		if (fabs(mx) > fabs(mn))
			scales[j] = mx/32767.0;
		else
			scales[j] = mn/32768.0;

		dptr[j] = (short *)mxCalloc(nk, sizeof(short));
		for(k=0; k<nk; k++)
		{
			dptr[j][k] = (short)rint(sptr[k]/scales[j]);
		}
	}

	scale = 0.0;
	for(j=0; j<nj; j++)
	{
		if (scales[j]>scale)
			scale = scales[j];
	}

	for(j=0; j<nj; j++)
	{
		for(k=0; k<nk; k++)
			dptr[j][k] = (short)rint(dptr[j][k]*(scales[j]/scale));
		write(fo, dptr[j], sizeof(short)*nk);
		mxFree((char *)(dptr[j]));
	}
	close(fo);

	mxFree((char *)scales);
	mxFree((char *)sptr);
	mxFree((char *)image);
	mxFree((char *)dptr);

	free_maps(maps);

	if (nlhs == 1)
	{
		plhs[0] = mxCreateFull(1,1, REAL);
		mxGetPr(plhs[0])[0] = scale;
	}
}
