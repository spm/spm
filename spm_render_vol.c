#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif

#include "cmex.h"
#include<math.h>

void surface(mat, zbuff, xcords, ycords, zcords, xdim1, ydim1, bytedat, xdim2, ydim2, zdim2, thresh)
double  mat[16];
double xcords[], ycords[], zcords[], zbuff[];
unsigned char bytedat[];
int xdim1, ydim1, xdim2, ydim2, zdim2, thresh;
{
	/*
	Project the coordinates of voxels which lie on the surface of the object
	onto the viewing plane.
	A surface is defined as a voxel above "thresh" adjascent to voxel(s) below "thresh"
	- assuming 6 nearest neighbours.
	Only the voxels closest to the viewing plane are stored.

	Projection performed by matrix "mat".
	The elements of "mat" are numbered:
		0  4  8 12
		1  5  9 13
		2  6 10 14
		3  7 11 15

	voxel (x1, y1, z1) is projected onto viewing plane at:
	((x1*mat[0] + y1*mat[4] + z1*mat[8] + mat[12]),
	(x1*mat[1] + y1*mat[5] + z1*mat[9] + mat[13]))
	Distance from the viewing plane (stored in "zbuff") is:
	(x1*mat[2] + y1*mat[6] + z1*mat[10] + mat[14])

	Note that elements 12-15 are unused, and that the projection is a parallel one.
	
	*/
	int z;
	double xs, ys;
	unsigned char *p;

	/* size of rectangle which is projected */
	xs = (0.01+sqrt(mat[0]*mat[0]+mat[4]*mat[4]+mat[8]*mat[8]));
	ys = (0.01+sqrt(mat[1]*mat[1]+mat[5]*mat[5]+mat[9]*mat[9]));

	p = bytedat;
	for(z=1; z<=zdim2; z++)
	{
		int y;
		double x2, y2, z2;
		x2 = mat[12] + z*mat[8] -1.0;
		y2 = mat[13] + z*mat[9] -1.0;
		z2 = mat[14] + z*mat[10];

		for(y=1; y<=ydim2; y++)
		{
			int x;
			double x3 = x2 + y*mat[4];
			double y3 = y2 + y*mat[5];
			double z3 = z2 + y*mat[6];

			for(x=1; x<=xdim2; x++)
			{
				/* check for surface voxel - this could be optimised greatly */
				if (*p > thresh && ((x<xdim2 && *(p+1) < thresh) || (x>1 && *(p-1) < thresh)
					|| (y<ydim2 && *(p+xdim2)<thresh) || (y>1 && *(p-xdim2)<thresh)
					|| (z<zdim2 && *(p+xdim2*ydim2)<thresh) || (z>1 && *(p-xdim2*ydim2)<thresh)))
				{
					int ix4, iy4, ixstart, ixend, iystart, iyend;
					double x4 = (x3 + x*mat[0]);
					double y4 = (y3 + x*mat[1]);
					double z4 = (z3 + x*mat[2]);
					ixstart = nint(x4);
					if (ixstart < 0) ixstart = 0;
					ixend = nint(x4+xs);
					if (ixend >= xdim1) ixend = xdim1-1;
					iystart = nint(y4);
					if (iystart < 0) iystart = 0;
					iyend = nint(y4+ys);
					if (iyend >= ydim1) iyend = ydim1-1;
					for(iy4 = iystart ; iy4<=iyend; iy4++)
						for(ix4 = ixstart ; ix4<=ixend; ix4++)
						{
							int off =  ix4+iy4*xdim1;
							/* is new voxel closer than original one */
							if (z4 < zbuff[off])
							{
								zbuff[off] = z4;
								xcords[off] = x;
								ycords[off] = y;
								zcords[off] = z;
							}
						}
				}
				p++;
			}
		}
	}
}

void render(xcords, ycords, zcords, xdim1, ydim1, out, bytedat, xdim2, ydim2, zdim2, light, thresh, nn)
double light[3];
double out[], xcords[], ycords[], zcords[];
unsigned char bytedat[];
int xdim1, ydim1, xdim2, ydim2, zdim2;
int thresh, nn;
{
	double	*filt, *filtp;
	int	dx, dy, dz, i;

	if ((filt = (double *)mxCalloc((nn*2+1)*(nn*2+1)*(nn*2+1), sizeof(double))) == (double *)0)
		return;

	/* create a 3D 'filter' which should produce the appropriate shading */
	filtp = filt;
	for(dz=-nn;dz<=nn;dz++)
			for(dy=-nn;dy<=nn;dy++)
				for(dx=-nn;dx<=nn;dx++)
					*(filtp++) = 
						(-light[0]*dx-light[1]*dy-light[2]*dz)/(nn*nn*nn*8);

	for(i=0; i<xdim1*ydim1; i++)
		if (xcords[i])
		{
			int x, y, z;
			double val;
			x = xcords[i]-1;
			y = ycords[i]-1;
			z = zcords[i]-1;
			filtp = filt;
			val = 0.0;
			for(dz=z-nn; dz<=z+nn; dz++)
				for(dy=y-nn; dy<=y+nn; dy++)
					for(dx=x-nn; dx<=x+nn; dx++)
					{
						if (dx>=0 && dx<xdim2 &&
							dy>=0 && dy<ydim2 && dz>=0 && dz<zdim2)
						{
							if (bytedat[dx+xdim2*(dy+ydim2*dz)]<thresh)
								val += *filtp;
						}
						filtp++;
					}
			if (val<0) val = 0.0;
			out[i] = val;
		}
	(void)mxFree((char *)filt);
}

void initdat(dat, size,val)
double dat[], val;
int size;
{
	int i;
	for(i=0; i<size; dat[i++] = val);
}

#include "volume.h"

#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	MAPTYPE *map;
	int m, n, xdim, ydim, zdim, k;
	double *mat;
	double *out, *odims, *params, *zbuff, *x, *y, *z;
	double light[3];

	if (nrhs != 4 || nlhs > 5) mexErrMsgTxt("Inappropriate usage.");

	map = get_map(prhs[0]);
	if (map->datatype != UNSIGNED_CHAR)
	{
		mexErrMsgTxt("Can only handle 8 bit volume.");
	}

	xdim = abs(nint(map->xdim));
	ydim = abs(nint(map->ydim));
	zdim = abs(nint(map->zdim));

	for(k=1; k<nrhs; k++)
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

	if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 2)
		mexErrMsgTxt("Output dimensions must have two elements.");
	odims = mxGetPr(prhs[2]);
	m = (int)odims[0];
	n = (int)odims[1];

	if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 2)
		mexErrMsgTxt("Parameters must have 2 elements.");
	params = mxGetPr(prhs[3]);

	plhs[0] = mxCreateFull(m,n, REAL);
	out = mxGetPr(plhs[0]);
	plhs[1] = mxCreateFull(m,n, REAL);
	zbuff = mxGetPr(plhs[1]);
	plhs[2] = mxCreateFull(m,n, REAL);
	x = mxGetPr(plhs[2]);
	plhs[3] = mxCreateFull(m,n, REAL);
	y = mxGetPr(plhs[3]);
	plhs[4] = mxCreateFull(m,n, REAL);
	z = mxGetPr(plhs[4]);

	initdat(zbuff,m*n,1024.0);
	surface(mat, zbuff, x, y, z, m, n, map->data, xdim, ydim, zdim, (int)params[0]);

	light[0] = mat[2];
	light[1] = mat[6];
	light[2] = mat[10];

	render(x, y, z, m, n, out, map->data, xdim, ydim, zdim, light, (int)params[0], (int)params[1]);
}
