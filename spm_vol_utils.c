#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif

#define PI 3.14159265358979
#define TINY 5e-2

#ifdef UNSIGNED_CHAR
#define RESAMPLE resample_uchar
#define SLICE slice_uchar
#define IMAGE_DTYPE unsigned char
#define RESAMPLE_0 resample_uchar_0
#define RESAMPLE_1 resample_uchar_1
#define RESAMPLE_SINC resample_uchar_sinc
#define SLICE_0 slice_uchar_0
#define SLICE_1 slice_uchar_1
#define SLICE_SINC slice_uchar_sinc

#elif SIGNED_SHORT
#define RESAMPLE resample_short
#define SLICE slice_short
#define IMAGE_DTYPE short int
#define RESAMPLE_0 resample_short_0
#define RESAMPLE_1 resample_short_1
#define RESAMPLE_SINC resample_short_sinc
#define SLICE_0 slice_short_0
#define SLICE_1 slice_short_1
#define SLICE_SINC slice_short_sinc

#elif SIGNED_INT
#define RESAMPLE resample_int
#define SLICE slice_int
#define IMAGE_DTYPE int
#define RESAMPLE_0 resample_int_0
#define RESAMPLE_1 resample_int_1
#define RESAMPLE_SINC resample_int_sinc
#define SLICE_0 slice_int_0
#define SLICE_1 slice_int_1
#define SLICE_SINC slice_int_sinc

#elif FLOAT
#define RESAMPLE resample_float
#define SLICE slice_float
#define IMAGE_DTYPE float
#define RESAMPLE_0 resample_float_0
#define RESAMPLE_1 resample_float_1
#define RESAMPLE_SINC resample_float_sinc
#define SLICE_0 slice_float_0
#define SLICE_1 slice_float_1
#define SLICE_SINC slice_float_sinc

#elif DOUBLE
#define RESAMPLE resample_double
#define SLICE slice_double
#define IMAGE_DTYPE double
#define RESAMPLE_0 resample_double_0
#define RESAMPLE_1 resample_double_1
#define RESAMPLE_SINC resample_double_sinc
#define SLICE_0 slice_double_0
#define SLICE_1 slice_double_1
#define SLICE_SINC slice_double_sinc

#else
#define RESAMPLE resample_uchar
#define SLICE slice_uchar
#define IMAGE_DTYPE unsigned char
#define RESAMPLE_0 resample_uchar_0
#define RESAMPLE_1 resample_uchar_1
#define RESAMPLE_SINC resample_uchar_sinc
#define SLICE_0 slice_uchar_0
#define SLICE_1 slice_uchar_1
#define SLICE_SINC slice_uchar_sinc
#endif

#include <math.h>

#ifndef NOLOOKUP
/* Generate a sinc lookup table with a Hanning filter envelope */
void make_lookup(coord,nn,dim, d1, table,ptpend)
double coord;
double table[], **ptpend;
int nn, dim, *d1;
{
	register int d2, d, fcoord;
	register double *tp, *tpend, dtmp;

	if (fabs(coord-rint(coord))<0.00001)
	{
		/* Close enough to use nearest neighbour */
		*d1=rint(coord);
		if (*d1<1 || *d1>dim) /* Pixel location outside image */
			*ptpend = table-1;
		else
		{
			table[0]=1.0;
			*ptpend = table;
		}
	}
	else
	{
		fcoord = floor(coord);
		*d1 = fcoord-nn;
		if (*d1<1) *d1=1;
		d2 = fcoord+nn;
		if (d2>dim) d2 = dim;

		*ptpend = tpend = table+(d2 - *d1);
		d = *d1, tp = table;
		while (tp <= tpend)
		{
			dtmp = PI*(coord-(d++));
			*(tp++) = sin(dtmp)/dtmp* 0.5*(1.0 + cos(dtmp/nn)) ;
		}
	}
}
#endif

/* Zero order hold resampling - nearest neighbour */
void RESAMPLE_0(m,vol,out,x,y,z,xdim,ydim,zdim,background, scale)
int m, xdim,ydim,zdim;
double out[], x[], y[], z[], background, scale;
IMAGE_DTYPE vol[];
{
	int i;
	vol -= (1 + xdim*(1 + ydim));
	for (i=0; i<m; i++)
	{
		int xcoord, ycoord, zcoord;
		xcoord = x[i]+0.5;
		ycoord = y[i]+0.5;
		zcoord = z[i]+0.5;
		if (xcoord>=1 && xcoord<=xdim && ycoord>=1 &&
			ycoord<=ydim && zcoord>=1 && zcoord<=zdim)
			out[i] = scale*vol[xcoord  + xdim*(ycoord  + ydim*(zcoord))];
		else out[i] = background;
	}
}


/* First order hold resampling - trilinear interpolation */
void RESAMPLE_1(m,vol,out,x,y,z,xdim,ydim,zdim,background, scale)
int m, xdim,ydim,zdim;
double out[], x[], y[], z[], background, scale;
IMAGE_DTYPE vol[];
{
	int i, dim1xdim2 = xdim*ydim;
	vol -= (1 + xdim*(1 + ydim));
	for (i=0; i<m; i++)
	{
		double xi,yi,zi;
		xi=x[i];
		yi=y[i];
		zi=z[i];

		if (zi>=1-TINY && zi<=zdim+TINY &&
			yi>=1-TINY && yi<=ydim+TINY
			&& xi>=1-TINY && xi<=xdim+TINY)
		{
			double k111,k112,k121,k122,k211,k212,k221,k222;
			double dx1, dx2, dy1, dy2, dz1, dz2;
			int off1, off2, offx, offy, offz, xcoord, ycoord, zcoord;

			xcoord = (int)floor(xi); dx1=xi-xcoord; dx2=1.0-dx1;
			ycoord = (int)floor(yi); dy1=yi-ycoord; dy2=1.0-dy1;
			zcoord = (int)floor(zi); dz1=zi-zcoord; dz2=1.0-dz1;

			xcoord = (xcoord < 1) ? ((offx=0),1) : ((offx=(xcoord>=xdim) ? 0 : 1        ),xcoord);
			ycoord = (ycoord < 1) ? ((offy=0),1) : ((offy=(ycoord>=ydim) ? 0 : xdim     ),ycoord);
			zcoord = (zcoord < 1) ? ((offz=0),1) : ((offz=(zcoord>=zdim) ? 0 : dim1xdim2),zcoord);

			off1 = xcoord  + xdim*(ycoord  + ydim*zcoord);
			k222 = vol[off1]; k122 = vol[off1+offx]; off2 = off1+offy;
			k212 = vol[off2]; k112 = vol[off2+offx]; off1+= offz;
			k221 = vol[off1]; k121 = vol[off1+offx]; off2 = off1+offy;
			k211 = vol[off2]; k111 = vol[off2+offx];

			/* resampled pixel value (trilinear interpolation) */
			out[i] = scale*(((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1)*dz2
				+ ((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1)*dz1);
		}
		else out[i] = background;

	}
}



/* Sinc resampling */
void RESAMPLE_SINC(m,vol,out,x,y,z,xdim,ydim,zdim, nn,background, scale)
int m, xdim,ydim,zdim, nn;
double out[], x[], y[], z[], background, scale;
IMAGE_DTYPE vol[];
{
	int i, dim1xdim2 = xdim*ydim;
	int dx1, dy1, dz1;
	static double tablex[255], tabley[255], tablez[255];

	vol -= (1 + xdim*(1 + ydim));
	for (i=0; i<m; i++)
	{
		if (z[i]>=1-TINY && z[i]<=zdim+TINY &&
			y[i]>=1-TINY && y[i]<=ydim+TINY &&
			x[i]>=1-TINY && x[i]<=xdim+TINY)
		{
			IMAGE_DTYPE *dp1;
			double dat=0.0, *tp1, *tp1end, *tp2end, *tp3end;

			make_lookup(x[i], nn, xdim, &dx1, tablex, &tp3end);
			make_lookup(y[i], nn, ydim, &dy1, tabley, &tp2end);
			make_lookup(z[i], nn, zdim, &dz1, tablez, &tp1end);

			tp1 = tablez;
			dy1 *= xdim;
			dp1 = vol + dim1xdim2*dz1;

			while(tp1 <= tp1end)
			{
				IMAGE_DTYPE *dp2 = dp1 + dy1;
				double dat2 = 0.0,
				*tp2 = tabley;
				while (tp2 <= tp2end)
				{
					register double dat3 = 0.0, *tp3 = tablex;
					register IMAGE_DTYPE *dp3 = dp2 + dx1;
					while(tp3 <= tp3end)
						dat3 += *(dp3++) * *(tp3++);
					dat2 += dat3 * *(tp2++);
					dp2 += xdim;
				}
				dat += dat2 * *(tp1++);
				dp1 += dim1xdim2;
			}
			out[i] = scale*dat;
		}
		else out[i] = background;
	}
}


/* Zero order hold resampling - nearest neighbour */
SLICE_0(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, background, scale)
double  mat[16], background, scale;
double image[];
IMAGE_DTYPE vol[];
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
			if (s3 == 0.0) return(-1);
			ix4 = rint((x3 += dx3)/s3);
			iy4 = rint((y3 += dy3)/s3);
			iz4 = rint((z3 += dz3)/s3);
			if (iz4>=1 && iz4<=zdim2 && iy4>=1 && iy4<=ydim2 && ix4>=1 && ix4<=xdim2)
				image[t] = scale*(double)vol[ix4 + xdim2*(iy4 + ydim2*iz4)];
			else image[t] = background;
			t++;
		}
	}
	return(0);
}

#define TINY 5e-2

/* First order hold resampling - trilinear interpolation */
SLICE_1(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, background, scale)
double  mat[16], background, scale;
double image[];
IMAGE_DTYPE vol[];
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
			double x4,y4,z4;
			s3 += ds3;
			if (s3 == 0.0) return(-1);
			x4=(x3 += dx3)/s3;
			y4=(y3 += dy3)/s3;
			z4=(z3 += dz3)/s3;

			if (z4>=1-TINY && z4<=zdim2+TINY &&
				y4>=1-TINY && y4<=ydim2+TINY
				&& x4>=1-TINY && x4<=xdim2+TINY)
			{
				double k111,k112,k121,k122,k211,k212,k221,k222;
				double dx1, dx2, dy1, dy2, dz1, dz2;
				int off1, off2, offx, offy, offz, ix4, iy4, iz4;

				ix4 = floor(x4); dx1=x4-ix4; dx2=1.0-dx1;
				iy4 = floor(y4); dy1=y4-iy4; dy2=1.0-dy1;
				iz4 = floor(z4); dz1=z4-iz4; dz2=1.0-dz1;

				ix4 = (ix4 < 1) ? ((offx=0),1) : ((offx=(ix4>=xdim2) ? 0 : 1        ),ix4);
				iy4 = (iy4 < 1) ? ((offy=0),1) : ((offy=(iy4>=ydim2) ? 0 : xdim2    ),iy4);
				iz4 = (iz4 < 1) ? ((offz=0),1) : ((offz=(iz4>=zdim2) ? 0 : dim1xdim2),iz4);

				off1 = ix4  + xdim2*(iy4  + ydim2*iz4);
				k222 = vol[off1]; k122 = vol[off1+offx]; off2 = off1+offy;
				k212 = vol[off2]; k112 = vol[off2+offx]; off1+= offz;
				k221 = vol[off1]; k121 = vol[off1+offx]; off2 = off1+offy;
				k211 = vol[off2]; k111 = vol[off2+offx];

				/* resampled pixel value (trilinear interpolation) */
				image[t] = scale*(((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1)*dz2
					  + ((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1)*dz1);
			}
			else image[t] = background; 
			t++;
		}
	}
	return(0);
}


/* Sinc resampling */
SLICE_SINC(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, nn, background, scale)
int ydim1,xdim1, xdim2,ydim2,zdim2, nn;
double image[], mat[], background, scale;
IMAGE_DTYPE vol[];
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
			double x4,y4,z4;
			s3 += ds3;
			if (s3 == 0.0) return(-1);
			x4=(x3 += dx3)/s3;
			y4=(y3 += dy3)/s3;
			z4=(z3 += dz3)/s3;

			if (z4>=1-TINY && z4<=zdim2+TINY &&
				y4>=1-TINY && y4<=ydim2+TINY &&
				x4>=1-TINY && x4<=xdim2+TINY)
			{
				IMAGE_DTYPE *dp1;
				double dat=0.0, *tp1, *tp1end, *tp2end, *tp3end;

				make_lookup(x4, nn, xdim2, &dx1, tablex, &tp3end);
				make_lookup(y4 , nn, ydim2, &dy1, tabley, &tp2end);
				make_lookup(z4, nn, zdim2, &dz1, tablez, &tp1end);

				tp1 = tablez;
				dy1 *= xdim2;
				dp1 = vol + dim1xdim2*dz1;

				while(tp1 <= tp1end)
				{
					IMAGE_DTYPE *dp2 = dp1 + dy1;
					double dat2 = 0.0,
					*tp2 = tabley;
					while (tp2 <= tp2end)
					{
						register double dat3 = 0.0, *tp3 = tablex;
						register IMAGE_DTYPE *dp3 = dp2 + dx1;
						while(tp3 <= tp3end)
							dat3 += *(dp3++) * *(tp3++);
						dat2 += dat3 * *(tp2++);
						dp2 += xdim2;
					}
					dat += dat2 * *(tp1++);
					dp1 += dim1xdim2;
				}
				*(image++) = scale*dat;
			}
			else *(image++) = background;
		}
	}
	return(0);
}


SLICE(mat, image, xdim1,ydim1, vol, xdim2,ydim2,zdim2, hold,background, scale)
int ydim1,xdim1, xdim2,ydim2,zdim2, hold;
double image[], mat[], background, scale;
IMAGE_DTYPE vol[];
{
	if (hold == 0)
		return(SLICE_0(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, background, scale));
	if (hold == 1)
		return(SLICE_1(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, background, scale));
	else
		return(SLICE_SINC(mat, image, xdim1, ydim1, vol, xdim2, ydim2, zdim2, hold, background, scale));
}


void RESAMPLE(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale)
int m, xdim,ydim,zdim, hold;
double out[], x[], y[], z[], background, scale;
IMAGE_DTYPE vol[];
{
	if (hold == 0)
		RESAMPLE_0(m,vol,out,x,y,z,xdim,ydim,zdim, background, scale);
	else if (hold == 1)
		RESAMPLE_1(m,vol,out,x,y,z,xdim,ydim,zdim, background, scale);
	else
		RESAMPLE_SINC(m,vol,out,x,y,z,xdim,ydim,zdim, hold, background, scale);
}
