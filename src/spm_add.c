/*
 * $Id: spm_add.c 1893 2008-07-08 15:05:40Z john $
 */
 
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
#include "mex.h"
#include "spm_mapping.h"
#define RINT(A) floor((A)+0.5)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MAPTYPE *maps = NULL, *get_maps();
	double *sptr = NULL, *scales = NULL, *image = NULL, scale;
	short **dptr = NULL;
	int ni, nj, nk, i, j, k;
	static double mat[] = {1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1};
	int mask0flag = 0, floatflag = 0;
	double NaN = mxGetNaN();
	mxArray *wplane_args[3];
	int maxval = 0, minval = 0;
	int dtype;

	if ((nrhs != 2 && nrhs != 3) || nlhs > 1)
		mexErrMsgTxt("Incorrect usage.");
 	if (nrhs == 3)
 	{
 		if (!mxIsChar(prhs[2]))
 			mexErrMsgTxt("Incorrect usage.");
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
 				/* if (buf[i] == 'f') floatflag = 1; */
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
		minval = -2147483647;
		floatflag = 0;
	}
	else
	{
		floatflag = 1;
	}

	nj = maps[0].dim[2];
	nk = maps[0].dim[0]*maps[0].dim[1];

	/* The compiler doesn't like this line - but I think it's OK */
	wplane_args[0] = (struct mxArray_tag *)prhs[1];
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
			/* compiler complains here about 5th arg - was "maps[i]" */
			slice(mat, image, maps[i].dim[0],maps[i].dim[1], &maps[i], 0, 0.0);
			if (mask0flag &&
				(maps[i].dtype == 2   || maps[i].dtype == 4    || maps[i].dtype == 8   ||
				 maps[i].dtype == 512 || maps[i].dtype == 1024 || maps[i].dtype == 2048))
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
				if (!floatflag && !mxIsFinite(sptr[k])) sptr[k] = 0.0;
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
				dptr[j][k] = (short)RINT(sptr[k]/scales[j]);
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
		scale = scale*32767/maxval; /* should really also use minval */

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


