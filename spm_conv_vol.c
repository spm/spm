#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif

#include <math.h>
#include "cmex.h"
#include "volume.h"

#include<stdio.h>


convxy(pl, xdim, ydim, filtx, filty, fxdim, fydim, xoff, yoff, out, buff, datatype)
unsigned char pl[];
int xdim, ydim, fxdim, fydim, xoff, yoff, datatype;
double out[], filtx[], filty[], buff[];
{
	double *rs;
	int x,y,k;
	rs = out;
	for(y=0; y<ydim; y++)
	{
		switch (datatype)
		{
			case UNSIGNED_CHAR:
				for(x=0; x<xdim; x++)
					buff[x] = ((unsigned char *)pl)[x+y*xdim];
				break;
			case SIGNED_SHORT:
				for(x=0; x<xdim; x++)
					buff[x] = ((short *)pl)[x+y*xdim];
				break;
			case SIGNED_INT:
				for(x=0; x<xdim; x++)
					buff[x] = ((int *)pl)[x+y*xdim];
				break;
			case FLOAT:
				for(x=0; x<xdim; x++)
					buff[x] = ((float *)pl)[x+y*xdim];
				break;
			case DOUBLE:
				for(x=0; x<xdim; x++)
					buff[x] = ((double *)pl)[x+y*xdim];
				break;
			default:
				mexErrMsgTxt("This should not happen.");
				break;
		}
		for(x=0; x<xdim; x++)
		{
			double sum1 = 0.0;
			int fstart, fend;
			fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
			fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

			for(k=fstart; k<fend; k++)
				sum1 += buff[x-xoff-k]*filtx[k];
			rs[x] = sum1;
		}
		rs += xdim;
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
		rs += xdim;
	}
}


convxyz(vol, xdim, ydim, zdim, filtx, filty, filtz, fxdim, fydim, fzdim, xoff, yoff, zoff, fp, datatype, ovol)
unsigned char vol[];
int xdim, ydim, zdim;
int fxdim, fydim, fzdim, xoff, yoff, zoff, datatype;
double filtx[], filty[], filtz[], *ovol;
FILE *fp;
{
	double *tmp, *buff, **sortedv, *obuf;
	int xy, z, k, fstart, fend, startz, endz;

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
			convxy(vol + (z*xdim*ydim*get_datasize(datatype)/8), xdim, ydim,
				filtx, filty, fxdim, fydim, xoff, yoff,
				tmp+((z%fzdim)*xdim*ydim), buff, datatype);
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
			  switch (datatype)
			    {
			    case UNSIGNED_CHAR:
			      for(xy=0; xy<xdim*ydim; xy++)
				{
				  if (obuf[xy] > 255.0) obuf[xy]=255.0;
				  if (obuf[xy] < 0.0) obuf[xy]=0.0;
				  ((unsigned char *)obuf)[xy] = obuf[xy] + 0.5;
				}
			      break;
			    case SIGNED_SHORT:
			      for(xy=0; xy<xdim*ydim; xy++)
				{
				  if (obuf[xy] > 32767.0) obuf[xy]=32767.0;
				  if (obuf[xy] < -32768.0) obuf[xy]=-32768.0;
				  ((short *)obuf)[xy] = floor(obuf[xy]+0.5);
				}
			      break;
			    case SIGNED_INT:
			      for(xy=0; xy<xdim*ydim; xy++)
				((int *)obuf)[xy] = floor(obuf[xy]+0.5);
			      break;
			    case FLOAT:
			      for(xy=0; xy<xdim*ydim; xy++)
				((float *)obuf)[xy] = obuf[xy];
			      break;
			    case DOUBLE:
			      /* nothing needed */
			      break;
			    default:
			      mexErrMsgTxt("This should not happen.");
			      break;
			    }

			  if (fwrite((char *)obuf,
				     get_datasize(datatype)/8, xdim*ydim, fp) != xdim*ydim)
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




#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
        MAPTYPE *map;
        int k, stlen;
	FILE *fp;
	double *offsets;
	char *str;
	int xdim,ydim,zdim,datatype,FileOffset=0;
	double *VolInfo, *oVol;
	unsigned char *Vol;

        if (nrhs < 6 || nrhs > 7 || nlhs > 0)
        {
                mexErrMsgTxt("Inappropriate usage.");
        }
	
	/* Extra parameter flags use of image in memory */
	if (nrhs == 7) {
	  if (!mxIsNumeric(prhs[6]) || mxIsComplex(prhs[6]) ||
	      !mxIsFull(prhs[6]) || !mxIsDouble(prhs[6])) 
	    mexErrMsgTxt("Bad volume info");
	  VolInfo = mxGetPr(prhs[0]);
	  xdim = (int)VolInfo[0];
	  ydim = (int)VolInfo[1];
	  zdim = (int)VolInfo[2];
	  if (xdim*ydim*zdim != mxGetM(prhs[6])*mxGetN(prhs[6]))
	    mexErrMsgTxt("Map dimension info doesn't match array");
	  Vol = (unsigned char *)mxGetPr(prhs[6]);
	  datatype = DOUBLE;
	}
	else {
	  map=get_map(prhs[0]);
	  xdim = (int)map->xdim;
	  ydim = (int)map->ydim;
	  zdim = (int)map->zdim;
	  Vol = map->data;
	  datatype = map->datatype;
	  FileOffset = map->off;
	}

	/* Non string-ness of destination img flags use of memory for output */
	if (mxIsString(prhs[1])) {
	  stlen = mxGetN(prhs[1]);
	  str = (char *)mxCalloc(stlen+1, sizeof(char));
	  mxGetString(prhs[1],str,stlen+1);
	  
	  if ((fp = fopen(str,"r+")) == (FILE *)0)
	    {
	      if ((fp = fopen(str,"w")) == (FILE *)0)
		{
		  mxFree(str);
		  mexErrMsgTxt("Cant open output file.");
		}
	    }
	  (void)fseek(fp, (long)FileOffset, 0);
	} else {
	  if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) ||
	      !mxIsFull(prhs[1]) || !mxIsDouble(prhs[1]) ||
	      mxGetM(prhs[1])*mxGetN(prhs[1]) != xdim*ydim*zdim)
	    mexErrMsgTxt("Bad output array");
	  fp = (FILE *)0;
	  oVol = (double *)mxGetPr(prhs[1]);
	}
	  

        for(k=2; k<=5; k++)
                if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
                        !mxIsFull(prhs[k]) || !mxIsDouble(prhs[k]))
                        mexErrMsgTxt("Functions must be numeric, real, full and double.");

	if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 3)
		mexErrMsgTxt("Offsets must have three values.");
	offsets = mxGetPr(prhs[5]);

	if (convxyz(Vol, abs(xdim), abs(ydim), abs(zdim),
		mxGetPr(prhs[2]), mxGetPr(prhs[3]), mxGetPr(prhs[4]),
		mxGetM(prhs[2])*mxGetN(prhs[2]),
		mxGetM(prhs[3])*mxGetN(prhs[3]),
		mxGetM(prhs[4])*mxGetN(prhs[4]),
		(int)floor(offsets[0]), (int)floor(offsets[1]), (int)floor(offsets[2]),
		fp, datatype, oVol) != 0)
	{
		(void)fclose(fp);
		mexErrMsgTxt("Error writing data.");

	}
	(void)fclose(fp);

}
