#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif

#include <math.h>
#include "spm_mapping.h"

static void convxy(out, xdim, ydim, filtx, filty, fxdim, fydim, xoff, yoff, buff)
int xdim, ydim, fxdim, fydim, xoff, yoff;
double out[], filtx[], filty[], buff[];
{
	int x,y,k;
	for(y=0; y<ydim; y++)
	{
		for(x=0; x<xdim; x++)
		{
			buff[x] = out[x+y*xdim];
			if (!mxIsFinite(buff[x]))
				buff[x] = 0.0;
		}
		for(x=0; x<xdim; x++)
		{
			double sum1 = 0.0;
			int fstart, fend;
			fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
			fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

			for(k=fstart; k<fend; k++)
				sum1 += buff[x-xoff-k]*filtx[k];
			out[x+y*xdim] = sum1;
		}
	}
	for(x=0; x<xdim; x++)
	{
		for(y=0; y<ydim; y++)
			buff[y] = out[x+y*xdim];

		for(y=0; y<ydim; y++)
		{
			double sum1 = 0.0;
			int fstart, fend;
			fstart = ((y-yoff >= ydim) ? y-ydim-yoff+1 : 0);
			fend = ((y-(yoff+fydim) < 0) ? y-yoff+1 : fydim);

			for(k=fstart; k<fend; k++)
				sum1 += buff[y-yoff-k]*filty[k];
			out[y*xdim+x] = sum1;
		}
	}
}


static int convxyz(vol, filtx, filty, filtz, fxdim, fydim, fzdim, xoff, yoff, zoff, oVol, wplane_args)
MAPTYPE *vol;
int fxdim, fydim, fzdim, xoff, yoff, zoff;
double filtx[], filty[], filtz[];
double *oVol;
mxArray *wplane_args[3];
{
	double *tmp, *buff, **sortedv, *obuf;
	int xy, z, k, fstart, fend, startz, endz;
	static double mat[] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
	int xdim, ydim, zdim;

	xdim = vol->dim[0];
	ydim = vol->dim[1];
	zdim = vol->dim[2];

	tmp = (double *)mxCalloc(xdim*ydim*fzdim,sizeof(double));
	buff = (double *)mxCalloc(((ydim>xdim) ? ydim : xdim),sizeof(double));
	sortedv = (double **)mxCalloc(fzdim, sizeof(double *));


	startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
	endz   = zdim+fzdim+zoff-1;

	for (z=startz; z<endz; z++)
	{
		double sum2 = 0.0;

		if (z >= 0 && z<zdim)
		{
			mat[14] = z+1.0;
			slice(mat, tmp+((z%fzdim)*xdim*ydim), vol->dim[0], vol->dim[1], vol, 0, 0);
			convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
				filtx, filty, fxdim, fydim, xoff, yoff, buff);
		}
		if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
		{
			fstart = ((z >= zdim) ? z-zdim+1 : 0);
			fend = ((z-fzdim < 0) ? z+1 : fzdim);

			for(k=0; k<fzdim; k++)
			{
				int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
				sortedv[k] = &(tmp[z1*xdim*ydim]);
			}

			for(k=fstart, sum2=0.0; k<fend; k++)
				sum2 += filtz[k];

			if (!oVol)
				obuf = mxGetPr(wplane_args[1]);
			else
				obuf = &oVol[(z-fzdim-zoff+1)*ydim*xdim];

			if (sum2)
			{
				for(xy=0; xy<xdim*ydim; xy++)
				{
					double sum1=0.0;
					for(k=fstart; k<fend; k++)
						sum1 += filtz[k]*sortedv[k][xy];

					obuf[xy] = sum1/sum2;
				}
			}
			else
				for(xy=0; xy<xdim*ydim; xy++)
					obuf[xy] = 0.0;

			if (!oVol)
			{
				mxGetPr(wplane_args[2])[0] = z-fzdim-zoff+2.0;
				mexCallMATLAB(0, NULL, 3, wplane_args, "spm_write_plane");
			}
		}
	}
	mxFree((char *)tmp);
	mxFree((char *)buff);
	mxFree((char *)sortedv);
	return(0);
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        MAPTYPE *map, *get_maps();
        int k;
	double *offsets, *oVol = NULL;
	mxArray *wplane_args[3];

	if (nrhs < 6 || nlhs > 0)
	{
		mexErrMsgTxt("Inappropriate usage.");
	}
	
	map=get_maps(prhs[0], &k);
	if (k!=1)
	{
		free_maps(map, k);
		mexErrMsgTxt("Too many images to smooth at once.");
	}

	if (!mxIsNumeric(prhs[1]))
	{
		/* The compiler doesn't like this line - but I think it's OK */
		wplane_args[0] = (struct mxArray_tag *)prhs[1];
		wplane_args[1] = mxCreateDoubleMatrix(map->dim[0],map->dim[1],mxREAL);
		wplane_args[2] = mxCreateDoubleMatrix(1,1,mxREAL);
		oVol = (double *)0;
	}
	else
	{
		if (   mxIsComplex(prhs[1]) ||
			mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) ||
			mxGetM(prhs[1])*mxGetN(prhs[1]) != map->dim[0]*map->dim[1]*map->dim[2])
		{
			free_maps(map, 1);
			mexErrMsgTxt("Bad output array");
		}
		oVol = (double *)mxGetPr(prhs[1]);
	}

        for(k=2; k<=5; k++)
	{
                if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
                        mxIsSparse(prhs[k]) || !mxIsDouble(prhs[k]))
		{
			free_maps(map, 1);
                        mexErrMsgTxt("Functions must be numeric, real, full and double.");
		}
	}

	if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 3)
	{
		free_maps(map, 1);
		mexErrMsgTxt("Offsets must have three values.");
	}
	offsets = mxGetPr(prhs[5]);

	if (convxyz(map,
		mxGetPr(prhs[2]), mxGetPr(prhs[3]), mxGetPr(prhs[4]),
		mxGetM(prhs[2])*mxGetN(prhs[2]),
		mxGetM(prhs[3])*mxGetN(prhs[3]),
		mxGetM(prhs[4])*mxGetN(prhs[4]),
		(int)floor(offsets[0]), (int)floor(offsets[1]), (int)floor(offsets[2]),
		oVol, wplane_args) != 0)
	{
		free_maps(map, 1);
		mexErrMsgTxt("Error writing data.");
	}
	free_maps(map, 1);
	if (!oVol)
	{
		mxDestroyArray(wplane_args[1]);
		mxDestroyArray(wplane_args[2]);
	}
}
