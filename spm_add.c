#ifndef lint
static char sccsid[]="%W% John Ashburner & Jean-Baptiste Poline%E%";
#endif
 
/*
% add a series of images - a compiled routine
% FORMAT s = spm_add(V,Q,flags)
% V     - Vector of mapped volumes (from spm_map or spm_vol).
% Q     - Filename for averaged image
% flags - Flags can be:
%               'f' - writes floating point output image.
%               'm' - masks the mean to zero or NaN wherever
%                     a zero occurs in the input images.
% s     - Scalefactor for output image.
%_______________________________________________________________________
%
% spm_add computes a sum of a set of image volumes to produce an 
% integral image that is written to a named file (Q).
%
% The image is written as signed short (16 bit) unless the `f' flag 
% is specified. 
%
% A mean can be effected by scaling the output image via it's
% scalefactor (see spm_mean for an example). A weighted sum can be
% effected by weighting the image scalefactors appropriately.
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
	int floatflag = 0;
	double NaN = 0.0/0.0; /* the only way to get a NaN that I know */
	mxArray *wplane_args[3];
	int maxval, minval, dtype;

	if ((nrhs != 2) || nlhs > 1)
		mexErrMsgTxt("Inappropriate usage.");

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

	dtype = get_dtype(prhs[1]);
	if (dtype > 256)
		dtype>>=8;

	if (dtype == 2)
	{
		maxval = 255;
		minval = 0;
		floatflag = 0;
	}
	else if (dtype == 4)
	{
		maxval = 32767;
		minval = -32768;
		floatflag = 0;
	}
	else if (dtype == 8)
	{
		maxval = 2147483647;
		minval = -2147483648;
		floatflag = 0;
	}
	else
	{
		floatflag = 1;
	}

	nj = maps[0].dim[2];
	nk = maps[0].dim[0]*maps[0].dim[1];

	/* The compiler doesn't like this line - but I think it's OK */
	wplane_args[0] = prhs[1];
	wplane_args[1] = mxCreateDoubleMatrix(maps[0].dim[0],maps[0].dim[1],mxREAL);
	wplane_args[2] = mxCreateDoubleMatrix(1,1,mxREAL);

	sptr   = mxGetPr(wplane_args[1]);
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
			if (maps[i].dtype == 2   || maps[i].dtype == 4    || maps[i].dtype == 8 ||
			    maps[i].dtype == 512 || maps[i].dtype == 1024 || maps[i].dtype == 2048)
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
			mxGetPr(wplane_args[2])[0] = j+1.0;
			mexCallMATLAB(0, NULL, 3, wplane_args, "spm_write_plane");
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
			if (scales[j] > scale)
				scale = scales[j];
		}
		scale = scale*32767/maxval;

		for(j=0; j<nj; j++)
		{
			for(k=0; k<nk; k++)
				sptr[k] = dptr[j][k]*(scales[j]/scale);

			mxGetPr(wplane_args[2])[0] = j+1.0;
			mexCallMATLAB(0, NULL, 3, wplane_args, "spm_write_plane");

			mxFree((char *)(dptr[j]));
		}

		mxFree((char *)scales);
		mxFree((char *)dptr);
	}
	else
	{
		scale = 1.0;
	}

	mxFree((char *)image);

	free_maps(maps, ni);

	if (nlhs == 1)
	{
		plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
		mxGetPr(plhs[0])[0] = scale;
	}
}


int get_dtype(const mxArray *ptr)
{
	mxArray *tmp;
	double *pr;

	if (!mxIsStruct(ptr))
	{
		mexErrMsgTxt("Not a structure.");
		return(NULL);
	}
	tmp=mxGetField(ptr,0,"dim");
	if (tmp == (mxArray *)0)
	{
		mexErrMsgTxt("Cant find dim.");
	}
	if (mxGetM(tmp)*mxGetN(tmp) != 4)
	{
		mexErrMsgTxt("Wrong sized dim.");
	}
	pr = mxGetPr(tmp);

	return((int)fabs(pr[3]));
}
