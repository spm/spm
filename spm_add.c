#ifndef lint
static char sccsid[]="%W% Jean-Baptiste Poline & John Ashburner %E%";
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
#include "spm_vol_utils.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MAPTYPE *maps, *get_maps();
	double *sptr, *scales, *image, scale;
	short **dptr;
	int ni, nj, nk, i, j, k, fo;
	static double mat[] = {1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1};
	char label[1024];
	int mask0flag = 0, floatflag = 0;
	double NaN = 0.0/0.0;

	if ((nrhs != 2 && nrhs != 3) || nlhs > 1)
		mexErrMsgTxt("Inappropriate usage.");

	if (nrhs == 3)
	{
		if (!mxIsChar(prhs[2]))
			mexErrMsgTxt("Inappropriate usage.");
		else
		{
			char *buf;
			int buflen;
			buflen = mxGetN(prhs[2])*mxGetM(prhs[2])+1;
			buf = mxCalloc(buflen,sizeof(char));
			if (mxGetString(prhs[2],buf,buflen))
			{
				mxFree(buf);
				mexErrMsgTxt("Cant get flags.");
			}
			for (i=0; i<buflen; i++)
			{
				if (buf[i] == 'm') mask0flag = 1;
				if (buf[i] == 'f') floatflag = 1;
			}
			mxFree(buf);
		}
	}

	maps = get_maps(prhs[0], &ni);


	for(i=1; i<ni; i++)
	{
		if (	maps[i].dim[0] != maps[0].dim[0] ||
			maps[i].dim[1] != maps[0].dim[1] ||
			maps[i].dim[2] != maps[0].dim[2])
			{
				free_maps(maps, ni);
				mexErrMsgTxt("Incompatible image dimensions.");
			}
	}

	mxGetString(prhs[1],label,mxGetN(prhs[1])+1);
	fo  = open(label, O_RDWR|O_CREAT, 0644);
	if (fo == -1)
		mexErrMsgTxt("Cant open output file.");

	nj = maps[0].dim[2];
	nk = maps[0].dim[0]*maps[0].dim[1];

	sptr   = (double *)mxCalloc(nk, sizeof(double));
	image  = (double *)mxCalloc(nk, sizeof(double));

	if (!floatflag)
	{
		scales = (double *)mxCalloc(nj, sizeof(double));
		dptr   = (short **)mxCalloc(nj, sizeof(double));
	}

	for(j=0; j<maps[0].dim[2]; j++)
	{
		double mx, mn;
		mat[14] = j+1.0;

		for(k=0; k<nk; k++)
			sptr[k] = 0.0;

		for(i=0; i<ni; i++)
		{
			slice(mat, image, maps[i].dim[0],maps[i].dim[1], maps[i], 0, 0.0);
			if (mask0flag)
			{
				for(k=0; k<nk; k++)
				{
					if (image[k] != 0)
						sptr[k] += image[k];
					else
						sptr[k] = NaN;
				}
			}
			else
			{
				for(k=0; k<nk; k++)
					sptr[k] += image[k];
			}
		}

		if (floatflag)
		{
			write(fo, sptr, sizeof(double)*nk);
		}
		else
		{
			/* Determine maximum and minimum */
			mx = -9e99;
			mn = 9e99;
			for(k=0; k<nk; k++)
			{
				if (!floatflag && !finite(sptr[k])) sptr[k] = 0.0;
				if (sptr[k]>mx) mx=sptr[k];
				if (sptr[k]<mn) mn=sptr[k];
			}

			if (mx > -mn)
				scales[j] = mx/32767.0;
			else
				scales[j] = -mn/32768.0;

			dptr[j] = (short *)mxCalloc(nk, sizeof(short));
			for(k=0; k<nk; k++)
			{
				dptr[j][k] = (short)rint(sptr[k]/scales[j]);
			}
		}
	}

	if (!floatflag)
	{

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

		mxFree((char *)scales);
		mxFree((char *)dptr);
	}
	else
	{
		scale = 1.0;
	}

	mxFree((char *)sptr);
	mxFree((char *)image);
	close(fo);

	free_maps(maps, ni);

	if (nlhs == 1)
	{
		plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
		mxGetPr(plhs[0])[0] = scale;
	}
}
