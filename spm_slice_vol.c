#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <math.h>
#include "cmex.h"
#include "volume.h"

/* Zero order hold resampling - nearest neighbour */
void slice_8bit_0(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2)
double  mat[16];
double image[];
unsigned char vol[];
int xdim1, ydim1, xdim2, ydim2, zdim2;
{
	double y, x2, y2, z2, s2, dx3=mat[0], dy3=mat[1], dz3=mat[2], ds3=mat[3];
	int t = 0;
	vol -= (1 + xdim2*(1 + ydim2));

	x2 = mat[12] + 0*mat[8];
	y2 = mat[13] + 0*mat[9];
	z2 = mat[14] + 0*mat[10];
	s2 = mat[15] + 0*mat[11];

	for(y=1; y<=ydim1; y++)
	{
		double x;
		double x3 = x2 + y*mat[4];
		double y3 = y2 + y*mat[5];
		double z3 = z2 + y*mat[6];
		double s3 = s2 + y*mat[7];

		for(x=1; x<=xdim1; x++)
		{
			int ix4, iy4, iz4;
			s3 += ds3;
			if (s3 == 0.0) mexErrMsgTxt("Fatal division by zero from transformation matrix.");
			ix4 = nint((x3 += dx3)/s3);
			iy4 = nint((y3 += dy3)/s3);
			iz4 = nint((z3 += dz3)/s3);
			if (iz4>=1 && iz4<=zdim2 && iy4>=1 && iy4<=ydim2 && ix4>=1 && ix4<=xdim2)
				image[t] = (double)vol[ix4 + xdim2*(iy4 + ydim2*iz4)];
			t++;
		}
	}
}


/* First order hold resampling - trilinear interpolation */
void slice_8bit_1(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2)
double  mat[16];
double image[];
unsigned char vol[];
int xdim1, ydim1, xdim2, ydim2, zdim2;
{
	double y, x2, y2, z2, s2, dx3=mat[0], dy3=mat[1], dz3=mat[2], ds3=mat[3];
	int t = 0, dim1xdim2 = xdim2*ydim2;
	vol -= (1 + xdim2*(1 + ydim2));

	x2 = mat[12] + 0*mat[8];
	y2 = mat[13] + 0*mat[9];
	z2 = mat[14] + 0*mat[10];
	s2 = mat[15] + 0*mat[11];

	for(y=1; y<=ydim1; y++)
	{
		double x;
		double x3 = x2 + y*mat[4];
		double y3 = y2 + y*mat[5];
		double z3 = z2 + y*mat[6];
		double s3 = s2 + y*mat[7];

		for(x=1; x<=xdim1; x++)
		{
			int ix4, iy4, iz4;
			double x4,y4,z4;
			s3 += ds3;
			if (s3 == 0.0) mexErrMsgTxt("Fatal division by zero from transformation matrix.");
			ix4 = (x4=(x3 += dx3)/s3);
			iy4 = (y4=(y3 += dy3)/s3);
			iz4 = (z4=(z3 += dz3)/s3);
			if (ix4>=1 && ix4<xdim2 && iy4>=1 && iy4<ydim2 && iz4>=1 && iz4<zdim2)
			{
				double k111,k112,k121,k122,k211,k212,k221,k222;
				double dx1, dx2, dy1, dy2, dz1, dz2;
				int off1, off2;
				dx1=x4-ix4; dx2=1.0-dx1;
				dy1=y4-iy4; dy2=1.0-dy1;
				dz1=z4-iz4; dz2=1.0-dz1;

				off1 = ix4  + xdim2*(iy4  + ydim2*iz4);
				k222 = vol[off1]; k122 = vol[off1+1]; off2 = off1+xdim2;
				k212 = vol[off2]; k112 = vol[off2+1]; off1+= dim1xdim2;
				k221 = vol[off1]; k121 = vol[off1+1]; off2 = off1+xdim2;
				k211 = vol[off2]; k111 = vol[off2+1];

				/* resampled pixel value (trilinear interpolation) */
				image[t] = nint(((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1)*dz2
					  + ((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1)*dz1);
			}
			t++;
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
	if (*d1<1)
		 *d1=1;
	d2 = fcoord+nn+1;
	if (d2>dim)
		d2 = dim;

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
void slice_8bit_sinc(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, nn)
int ydim1,xdim1, xdim2,ydim2,zdim2, nn;
double image[], mat[];
unsigned char vol[];
{
	int dim1xdim2 = xdim2*ydim2;
	int dx1, dy1, dz1;
	static double tablex[255], tabley[255], tablez[255];
	double y, x2, y2, z2, s2, dx3=mat[0], dy3=mat[1], dz3=mat[2], ds3=mat[3];

	vol -= (1 + xdim2*(1 + ydim2));

	x2 = mat[12] + 0*mat[8];
	y2 = mat[13] + 0*mat[9];
	z2 = mat[14] + 0*mat[10];
	s2 = mat[15] + 0*mat[11];

	for(y=1; y<=ydim1; y++)
	{
		double x;
		double x3 = x2 + y*mat[4];
		double y3 = y2 + y*mat[5];
		double z3 = z2 + y*mat[6];
		double s3 = s2 + y*mat[7];

		for(x=1; x<=xdim1; x++)
		{
			unsigned char *dp1;
			double dat=0.0, *tp1, *tp1end, *tp2end, *tp3end;

			s3 += ds3;
			if (s3 == 0.0) mexErrMsgTxt("Fatal division by zero from transformation matrix.");
			make_lookup((x3 += dx3)/s3, nn, xdim2, &dx1, tablex, &tp3end);
			make_lookup((y3 += dy3)/s3, nn, ydim2, &dy1, tabley, &tp2end);
			make_lookup((z3 += dz3)/s3, nn, zdim2, &dz1, tablez, &tp1end);

			tp1 = tablez;
			dy1 *= xdim2;
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
					dp2 += xdim2;
				}
				dat += dat2 * *(tp1++);
				dp1 += dim1xdim2;
			}
			*(image++) = dat;
		}
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
	MAPTYPE *map;
	int m,n, k, hold, xdim, ydim, zdim, datasize=8;
	double *mat, *ptr, scalefactor, *img;

	if (nrhs != 4 || nlhs > 1)
	{
		mexErrMsgTxt("Inappropriate usage.");
	}

	map = get_map(prhs[0]);

	if (map->datatype != UNSIGNED_CHAR)
	{
		mexErrMsgTxt("Can only handle 8 bit volume.");
	}
	datasize = 8;
	xdim = abs(nint(map->xdim));
	ydim = abs(nint(map->ydim));
	zdim = abs(nint(map->zdim));

	for(k=1; k<=3; k++)
	{
		if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
			!mxIsFull(prhs[k]) || !mxIsDouble(prhs[k]))
		{
			mexErrMsgTxt("Arguments must be numeric, real, full and double.");
		}
	}

	/* get transformation matrix */
	if (mxGetM(prhs[1]) != 4 && mxGetN(prhs[1]) != 4)
	{
		mexErrMsgTxt("Transformation matrix must be 4 x 4.");
	}
	mat = mxGetPr(prhs[1]);

	/* get output dimensions */
	if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != 2)
	{
		mexErrMsgTxt("Output dimensions must have two elements.");
	}
	ptr = mxGetPr(prhs[2]);
	m = abs(nint(ptr[0]));
	n = abs(nint(ptr[1]));
	plhs[0] = mxCreateFull(m,n,REAL);
	img = mxGetPr(plhs[0]);

	if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 1)
	{
		mexErrMsgTxt("Hold argument must have one element.");
	}
	hold = nint(*(mxGetPr(prhs[3])));

	if (hold == 0)
		slice_8bit_0(mat, img, m, n, (unsigned char *)(map->map), xdim, ydim, zdim);
	else if (hold == 1)
		slice_8bit_1(mat, img, m, n, (unsigned char *)(map->map), xdim, ydim, zdim);
	else if (hold >= 3 && hold <= 127)
		slice_8bit_sinc(mat, img, m, n, (unsigned char *)(map->map), xdim, ydim, zdim, hold);
	else
		mexErrMsgTxt("Bad hold value.");

	/* final rescale */
	if (map->scalefactor != 1.0 && map->scalefactor != 0)
		for (k=0; k<m*n; k++)
			img[k] *= map->scalefactor;
}
