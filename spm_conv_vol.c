#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif

#include <math.h>
#include "mex.h"
#include "spm_vol_utils.h"
#include "spm_map.h"

#include<stdio.h>


static void convxy(out, xdim, ydim, filtx, filty, fxdim, fydim, xoff, yoff, buff)
int xdim, ydim, fxdim, fydim, xoff, yoff;
double out[], filtx[], filty[], buff[];
{
	int x,y,k;
	for(y=0; y<ydim; y++)
	{
		for(x=0; x<xdim; x++)
			buff[x] = out[x+y*xdim];
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


static int convxyz(vol, filtx, filty, filtz, fxdim, fydim, fzdim, xoff, yoff, zoff, fp, ovol)
MAPTYPE *vol;
int fxdim, fydim, fzdim, xoff, yoff, zoff;
double filtx[], filty[], filtz[], *ovol;
FILE *fp;
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
	obuf = (double *)mxCalloc(xdim*ydim,sizeof(double));
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

			if (fp != (FILE *)0) {
			  switch (vol->dtype)
			    {
			    case UNSIGNED_CHAR:
			      for(xy=0; xy<xdim*ydim; xy++)
				{
				  obuf[xy] = (obuf[xy]-vol->offset)/vol->scale;
				  if (obuf[xy] > 255.0) obuf[xy]=255.0;
				  if (obuf[xy] < 0.0) obuf[xy]=0.0;
				  ((unsigned char *)obuf)[xy] = obuf[xy] + 0.5;
				}
			      break;
			    case SIGNED_SHORT:
			      for(xy=0; xy<xdim*ydim; xy++)
				{
				  obuf[xy] = (obuf[xy]-vol->offset)/vol->scale;
				  if (obuf[xy] > 32767.0) obuf[xy]=32767.0;
				  if (obuf[xy] < -32768.0) obuf[xy]=-32768.0;
				  ((short *)obuf)[xy] = floor(obuf[xy]+0.5);
				}
			      break;
			    case SIGNED_INT:
			      for(xy=0; xy<xdim*ydim; xy++)
			      {
				  obuf[xy] = (obuf[xy]-vol->offset)/vol->scale;
				  ((int *)obuf)[xy] = floor(obuf[xy]+0.5);
			      }
			      break;
			    case FLOAT:
			      for(xy=0; xy<xdim*ydim; xy++)
			      {
				  obuf[xy] = (obuf[xy]-vol->offset)/vol->scale;
				((float *)obuf)[xy] = obuf[xy];
			      }
			      break;
			    case DOUBLE:
			      for(xy=0; xy<xdim*ydim; xy++)
			      {
				  obuf[xy] = (obuf[xy]-vol->offset)/vol->scale;
			          /* nothing needed */
			      }
			      break;
			    default:
			      mexErrMsgTxt("This should not happen.");
			      break;
			    }

			  if (fwrite((char *)obuf,
				     get_datasize(vol->dtype)/8, xdim*ydim, fp) != xdim*ydim)
			    {
			      mxFree((char *)tmp);
			      mxFree((char *)buff);
			      mxFree((char *)obuf);
			      mxFree((char *)sortedv);
			      return(-1);
			    }

			} else {
			  for(xy=0; xy<xdim*ydim; xy++)
			    *(ovol++) = obuf[xy];
			}
		}
	}
	mxFree((char *)tmp);
	mxFree((char *)buff);
	mxFree((char *)obuf);
	mxFree((char *)sortedv);
	return(0);
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        MAPTYPE *map, *get_maps(), vol;
        int k, stlen;
	FILE *fp;
	double *offsets;
	char *str;
	int FileOffset=0;
	double *VolInfo, *oVol;

	if (nrhs < 6 || nrhs > 7 || nlhs > 0)
	{
		mexErrMsgTxt("Inappropriate usage.");
	}
	
	/* Extra parameter flags use of image in memory */
	if (nrhs == 7)
	{
		if (!mxIsNumeric(prhs[6]) || mxIsComplex(prhs[6]) ||
			mxIsSparse(prhs[6]) || !mxIsDouble(prhs[6])) 
			mexErrMsgTxt("Bad volume info");
		VolInfo = mxGetPr(prhs[0]);
		vol.dim[0] = (int)VolInfo[0];
		vol.dim[1] = (int)VolInfo[1];
		vol.dim[2] = (int)VolInfo[2];
		vol.scale  = 1.0;
		vol.offset = 0.0;
		vol.dtype  = DOUBLE;
		if (vol.dim[0]*vol.dim[1]*vol.dim[2] != mxGetM(prhs[6])*mxGetN(prhs[6]))
			mexErrMsgTxt("Map dimension info doesn't match array");
		vol.data = (unsigned char *)mxGetPr(prhs[6]);
		map = &vol;
	}
	else
	{
		MAP *lmap;
		map=get_maps(prhs[0], &k);
		if (k!=1)
		{
			free_maps(map);
			mexErrMsgTxt("Too many images to smooth at once.");
		}
		lmap=(MAP *)mxGetPr(prhs[0]);
		FileOffset = lmap->off;
	}

	/* Non string-ness of destination img flags use of memory for output */
	if (!mxIsNumeric(prhs[1]))
	{
		stlen = mxGetN(prhs[1]);
		str = (char *)mxCalloc(stlen+1, sizeof(char));
		if (mxGetString(prhs[1],str,stlen+1))
		{
			mxFree(str);
			mexErrMsgTxt("Could not convert string data.");
		}
		if ((fp = fopen(str,"r+")) == (FILE *)0)
		{
			if ((fp = fopen(str,"w")) == (FILE *)0)
			{
				mxFree(str);
				mexErrMsgTxt("Cant open output file.");
			}
		}
		(void)fseek(fp, (long)FileOffset, 0);
	}
	else
	{
		if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) ||
			mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) ||
			mxGetM(prhs[1])*mxGetN(prhs[1]) != map->dim[0]*map->dim[1]*map->dim[2])
			mexErrMsgTxt("Bad output array");
		fp = (FILE *)0;
		oVol = (double *)mxGetPr(prhs[1]);
	}
	  

        for(k=2; k<=5; k++)
                if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
                        mxIsSparse(prhs[k]) || !mxIsDouble(prhs[k]))
		{
                        mexErrMsgTxt("Functions must be numeric, real, full and double.");
		}

	if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 3)
	{
		mexErrMsgTxt("Offsets must have three values.");
	}
	offsets = mxGetPr(prhs[5]);

	if (convxyz(map,
		mxGetPr(prhs[2]), mxGetPr(prhs[3]), mxGetPr(prhs[4]),
		mxGetM(prhs[2])*mxGetN(prhs[2]),
		mxGetM(prhs[3])*mxGetN(prhs[3]),
		mxGetM(prhs[4])*mxGetN(prhs[4]),
		(int)floor(offsets[0]), (int)floor(offsets[1]), (int)floor(offsets[2]),
		fp, oVol) != 0)
	{
		(void)fclose(fp);
		mexErrMsgTxt("Error writing data.");
	}
	(void)fclose(fp);
}
