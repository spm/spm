#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <math.h>
#include "cmex.h"
#include "dbh.h"
#define MAGIC 170372
typedef struct mayomaptype
{
	int magic;
	struct dsr hdr;
	caddr_t map;
	size_t len;
	int prot;
	int flags;
}	MAYOMAPTYPE;


/* Zero order hold resampling - nearest neighbour */
void resample_8bit_0(m,vol,out,x,y,z,xdim,ydim,zdim)
int m, xdim,ydim,zdim;
double out[], x[], y[], z[];
unsigned char vol[];
{
	int i, dim1xdim2 = xdim*ydim;
	for (i=0; i<m; i++)
	{
		int xcoord, ycoord, zcoord;
		xcoord = nint(x[i]);
		ycoord = nint(y[i]);
		zcoord = nint(z[i]);
		if (xcoord>=1 && xcoord<=xdim && ycoord>=1 &&
			ycoord<=ydim && zcoord>=1 && zcoord<=zdim)
			out[i] = vol[xcoord-1  + xdim*(ycoord-1  + ydim*(zcoord-1))];
	}
}


/* First order hold resampling - trilinear interpolation */
void resample_8bit_1(m,vol,out,x,y,z,xdim,ydim,zdim)
int m, xdim,ydim,zdim;
double out[], x[], y[], z[];
unsigned char vol[];
{
	int i, dim1xdim2 = xdim*ydim;
	for (i=0; i<m; i++)
	{
		double xi,yi,zi,dx1,dx2,dy1,dy2,dz1,dz2,
			k111,k112,k121,k122,k211,k212,k221,k222;
		int xcoord, ycoord, zcoord, off1, off2;
		xi=x[i];
		yi=y[i];
		zi=z[i];

		xcoord = (int)floor(xi);
		ycoord = (int)floor(yi);
		zcoord = (int)floor(zi);
		if (xcoord>=1 && xcoord<xdim && ycoord>=1 &&
			ycoord<ydim && zcoord>=1 && zcoord<zdim)
		{
			dx1=xi-xcoord; dx2=1.0-dx1;
			dy1=yi-ycoord; dy2=1.0-dy1;
			dz1=zi-zcoord; dz2=1.0-dz1;

			off1 = xcoord-1  + xdim*(ycoord-1  + ydim*(zcoord-1));
			k222 = vol[off1]; k122 = vol[off1+1]; off2 = off1+xdim;
			k212 = vol[off2]; k112 = vol[off2+1]; off1+= dim1xdim2;
			k221 = vol[off1]; k121 = vol[off1+1]; off2 = off1+xdim;
			k211 = vol[off2]; k111 = vol[off2+1];

			/* resampled pixel value (trilinear interpolation) */
			out[i] = ((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1)*dz2
				+ ((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1)*dz1;
		}
	}
}

#define PI 3.14159265358979

void make_lookup(coord,nn,dim, d1, table,ptpend)
double coord;
double table[], **ptpend;
int nn, dim, *d1;
{
	register int d2, d, fcoord;
	register double *tp, *tpend, dtmp;
	fcoord = floor(coord);
	*d1 = fcoord-nn;
	if (*d1<1) *d1=1;
	d2 = fcoord+nn+1;
	if (d2>dim) d2 = dim;

	*ptpend = tpend = table+(d2 - *d1);
	d = *d1, tp = table;
	while (tp< tpend)
	{
		dtmp = PI*(coord-(d++));
		if (dtmp != 0)
			*(tp++) = sin(dtmp)/dtmp* 0.5*(1.0 + cos(dtmp/nn)) ;
		else
			*(tp++) = 1.0;
	}
}

/* Sinc resampling */
void resample_8bit_sinc(m,vol,out,x,y,z,xdim,ydim,zdim, nn)
int m, xdim,ydim,zdim, nn;
double out[], x[], y[], z[];
unsigned char vol[];
{
	int i, dim1xdim2 = xdim*ydim;
	int dx1, dy1, dz1;
	static double tablex[255], tabley[255], tablez[255];

	vol -= (1 + xdim*(1 + ydim));
	for (i=0; i<m; i++)
	{
		unsigned char *dp1;
		double dat=0.0, *tp1, *tp1end, *tp2end, *tp3end;

		make_lookup(x[i], nn, xdim, &dx1, tablex, &tp3end);
		make_lookup(y[i], nn, ydim, &dy1, tabley, &tp2end);
		make_lookup(z[i], nn, zdim, &dz1, tablez, &tp1end);

		tp1 = tablez;
		dy1 *= xdim;
		dp1 = vol + dim1xdim2*dz1;

		while(tp1<tp1end)
		{
			unsigned char *dp2 = dp1 + dy1;
			double dat2 = 0.0,
			*tp2 = tabley;
			while (tp2 < tp2end)
			{
				register double dat3 = 0.0, *tp3 = tablex;
				register unsigned char *dp3 = dp2 + dx1;
				while(tp3 < tp3end)
					dat3 += *(dp3++) * *(tp3++);
				dat2 += dat3 * *(tp2++);
				dp2 += xdim;
			}
			dat += dat2 * *(tp1++);
			dp1 += dim1xdim2;
		}
		out[i] = dat;
	}
}


#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	MAYOMAPTYPE *mayomap;
	struct image_dimension *dime;
	int m,n, k, hold;
	if (nrhs != 5 || nlhs > 1)
	{
		mexErrMsgTxt("Inappropriate usage.");
	}

	if ((!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ||
		!mxIsFull(prhs[0]) || !mxIsDouble(prhs[0]) ||
		mxGetN(prhs[0]) != 1 ||
		mxGetM(prhs[0]) != (sizeof(MAYOMAPTYPE)+sizeof(double)-1)/sizeof(double)))
	{
		mexErrMsgTxt("Bad image handle.");
	}
	mayomap = (MAYOMAPTYPE *)mxGetPr(prhs[0]);
	if (mayomap->magic != MAGIC)
	{
		mexErrMsgTxt("Bad magic number in image handle.");
	}

	dime = &(mayomap->hdr.dime);
	if (dime->bitpix != 8)
	{
		mexErrMsgTxt("Can only handle 8 bit volume.");
	}

	for(k=1; k<=3; k++)
	{
		if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
			!mxIsFull(prhs[k]) || !mxIsDouble(prhs[k]))
		{
			mexErrMsgTxt("Coordinates must be numeric, real, full and double.");
		}
	}

	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if (mxGetM(prhs[2]) != m || mxGetN(prhs[2]) != n ||
		mxGetM(prhs[2]) != m || mxGetN(prhs[2]) != n)
	{
		mexErrMsgTxt("Coordinates must have compatible dimensions.");
	}
	if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) ||
		!mxIsFull(prhs[4]) || !mxIsDouble(prhs[4]) ||
		mxGetM(prhs[4])*mxGetN(prhs[4]) != 1)
	{
		mexErrMsgTxt("Bad hold argument.");
	}
	hold = nint(*(mxGetPr(prhs[4])));
	plhs[0] = mxCreateFull(m,n,REAL);

	if (hold == 0)
		resample_8bit_0(m*n,(unsigned char *)(mayomap->map),mxGetPr(plhs[0]),
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]),
			dime->dim[1], dime->dim[2], dime->dim[3]);
	else if (hold == 1)
		resample_8bit_1(m*n,(unsigned char *)(mayomap->map),mxGetPr(plhs[0]),
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]),
			dime->dim[1], dime->dim[2], dime->dim[3]);
	else if (hold >= 3 && hold <= 127)
		resample_8bit_sinc(m*n,(unsigned char *)(mayomap->map),mxGetPr(plhs[0]),
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]),
			dime->dim[1], dime->dim[2], dime->dim[3], hold);
	else
		mexErrMsgTxt("Bad hold value.");
}
