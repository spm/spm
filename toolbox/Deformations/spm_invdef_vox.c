#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include "mex.h"
#include<math.h>
#include<stdio.h>
#define MAXV 2048
#define REAL float


static void invertX(REAL X[4][3], REAL IX[4][4], REAL dIX[4][3], REAL *dt)
{
	REAL id;
	*dt = X[0][0]*(X[3][1]*(X[1][2]-X[2][2])+X[1][1]*(X[2][2]-X[3][2])+X[2][1]*(X[3][2]-X[1][2]))+
	      X[1][0]*(X[3][2]*(X[0][1]-X[2][1])+X[0][2]*(X[2][1]-X[3][1])+X[2][2]*(X[3][1]-X[0][1]))+
	      X[2][0]*(X[0][1]*(X[1][2]-X[3][2])+X[3][1]*(X[0][2]-X[1][2])+X[1][1]*(X[3][2]-X[0][2]))+
	      X[3][0]*(X[1][2]*(X[2][1]-X[0][1])+X[0][2]*(X[1][1]-X[2][1])+X[2][2]*(X[0][1]-X[1][1]));
	id = 1/(*dt);
	IX[0][0] = id*(X[1][1]*(X[2][2]-X[3][2])+X[2][1]*(X[3][2]-X[1][2])+X[3][1]*(X[1][2]-X[2][2]));
	IX[0][1] = id*(X[0][1]*(X[3][2]-X[2][2])+X[2][1]*(X[0][2]-X[3][2])+X[3][1]*(X[2][2]-X[0][2]));
	IX[0][2] = id*(X[0][1]*(X[1][2]-X[3][2])+X[1][1]*(X[3][2]-X[0][2])+X[3][1]*(X[0][2]-X[1][2]));
	IX[0][3] = id*(X[0][1]*(X[2][2]-X[1][2])+X[1][1]*(X[0][2]-X[2][2])+X[2][1]*(X[1][2]-X[0][2]));
	IX[1][0] = id*(X[1][0]*(X[3][2]-X[2][2])+X[2][0]*(X[1][2]-X[3][2])+X[3][0]*(X[2][2]-X[1][2]));
	IX[1][1] = id*(X[0][0]*(X[2][2]-X[3][2])+X[2][0]*(X[3][2]-X[0][2])+X[3][0]*(X[0][2]-X[2][2]));
	IX[1][2] = id*(X[0][0]*(X[3][2]-X[1][2])+X[1][0]*(X[0][2]-X[3][2])+X[3][0]*(X[1][2]-X[0][2]));
	IX[1][3] = id*(X[0][0]*(X[1][2]-X[2][2])+X[1][0]*(X[2][2]-X[0][2])+X[2][0]*(X[0][2]-X[1][2]));
	IX[2][0] = id*(X[1][0]*(X[2][1]-X[3][1])+X[2][0]*(X[3][1]-X[1][1])+X[3][0]*(X[1][1]-X[2][1]));
	IX[2][1] = id*(X[0][0]*(X[3][1]-X[2][1])+X[2][0]*(X[0][1]-X[3][1])+X[3][0]*(X[2][1]-X[0][1]));
	IX[2][2] = id*(X[0][0]*(X[1][1]-X[3][1])+X[1][0]*(X[3][1]-X[0][1])+X[3][0]*(X[0][1]-X[1][1]));
	IX[2][3] = id*(X[0][0]*(X[2][1]-X[1][1])+X[1][0]*(X[0][1]-X[2][1])+X[2][0]*(X[1][1]-X[0][1]));
	IX[3][0] = id*(X[1][0]*(X[2][2]*X[3][1]-X[3][2]*X[2][1])+
		       X[2][0]*(X[3][2]*X[1][1]-X[1][2]*X[3][1])+
		       X[3][0]*(X[1][2]*X[2][1]-X[2][2]*X[1][1]));
	IX[3][1] = id*(X[0][0]*(X[3][2]*X[2][1]-X[2][2]*X[3][1])+
		       X[2][0]*(X[0][2]*X[3][1]-X[3][2]*X[0][1])+
		       X[3][0]*(X[2][2]*X[0][1]-X[0][2]*X[2][1]));
	IX[3][2] = id*(X[0][0]*(X[1][2]*X[3][1]-X[3][2]*X[1][1])+
		       X[1][0]*(X[3][2]*X[0][1]-X[0][2]*X[3][1])+
		       X[3][0]*(X[0][2]*X[1][1]-X[1][2]*X[0][1]));
	IX[3][3] = id*(X[0][0]*(X[2][2]*X[1][1]-X[1][2]*X[2][1])+
		       X[1][0]*(X[0][2]*X[2][1]-X[2][2]*X[0][1])+
		       X[2][0]*(X[1][2]*X[0][1]-X[0][2]*X[1][1]));

	dIX[0][0] = -( X[1][1]*X[2][2]-X[1][1]*X[3][2]-X[1][2]*X[2][1]+X[1][2]*X[3][1]+X[2][1]*X[3][2]-X[3][1]*X[2][2])*id;
	dIX[0][1] = -(-X[1][0]*X[2][2]+X[1][0]*X[3][2]+X[1][2]*X[2][0]-X[1][2]*X[3][0]-X[2][0]*X[3][2]+X[3][0]*X[2][2])*id;
	dIX[0][2] = -( X[1][0]*X[2][1]-X[1][0]*X[3][1]-X[1][1]*X[2][0]+X[1][1]*X[3][0]+X[2][0]*X[3][1]-X[3][0]*X[2][1])*id;

	dIX[1][0] =  ( X[0][1]*X[2][2]-X[0][1]*X[3][2]-X[0][2]*X[2][1]+X[0][2]*X[3][1]+X[2][1]*X[3][2]-X[3][1]*X[2][2])*id;
	dIX[1][1] =  ( X[3][0]*X[2][2]-X[2][0]*X[3][2]-X[0][0]*X[2][2]+X[0][0]*X[3][2]-X[0][2]*X[3][0]+X[0][2]*X[2][0])*id;
	dIX[1][2] =  ( X[0][0]*X[2][1]-X[0][0]*X[3][1]-X[0][1]*X[2][0]+X[0][1]*X[3][0]+X[2][0]*X[3][1]-X[3][0]*X[2][1])*id;

	dIX[2][0] = -( X[0][1]*X[1][2]-X[0][1]*X[3][2]-X[0][2]*X[1][1]+X[0][2]*X[3][1]+X[1][1]*X[3][2]-X[1][2]*X[3][1])*id;
	dIX[2][1] = -(-X[0][0]*X[1][2]+X[0][0]*X[3][2]+X[0][2]*X[1][0]-X[0][2]*X[3][0]-X[1][0]*X[3][2]+X[1][2]*X[3][0])*id;
	dIX[2][2] = -( X[0][0]*X[1][1]-X[0][0]*X[3][1]-X[0][1]*X[1][0]+X[0][1]*X[3][0]+X[1][0]*X[3][1]-X[1][1]*X[3][0])*id;

	dIX[3][0] =  ( X[0][1]*X[1][2]-X[0][1]*X[2][2]-X[0][2]*X[1][1]+X[0][2]*X[2][1]+X[1][1]*X[2][2]-X[1][2]*X[2][1])*id;
	dIX[3][1] =  ( X[1][2]*X[2][0]-X[1][0]*X[2][2]+X[0][0]*X[2][2]-X[0][0]*X[1][2]-X[0][2]*X[2][0]+X[0][2]*X[1][0])*id;
	dIX[3][2] =  ( X[0][0]*X[1][1]-X[0][0]*X[2][1]-X[0][1]*X[1][0]+X[0][1]*X[2][0]+X[1][0]*X[2][1]-X[1][1]*X[2][0])*id;
}

static void getM(REAL Y[4][3], REAL IX[4][4], REAL dIX[4][3], REAL M[4][3], int i, int j, int k, REAL *pix30)
{
	REAL ix30, ix31, ix32, ix33;

	ix30 = IX[3][0] + i*dIX[0][0] + j*dIX[0][1] + k*dIX[0][2];
	ix31 = IX[3][1] + i*dIX[1][0] + j*dIX[1][1] + k*dIX[1][2];
	ix32 = IX[3][2] + i*dIX[2][0] + j*dIX[2][1] + k*dIX[2][2];
	ix33 = IX[3][3] + i*dIX[3][0] + j*dIX[3][1] + k*dIX[3][2];

	M[0][0] = Y[0][0]*IX[0][0] + Y[1][0]*IX[0][1] + Y[2][0]*IX[0][2] + Y[3][0]*IX[0][3];
	M[0][1] = Y[0][1]*IX[0][0] + Y[1][1]*IX[0][1] + Y[2][1]*IX[0][2] + Y[3][1]*IX[0][3];
	M[0][2] = Y[0][2]*IX[0][0] + Y[1][2]*IX[0][1] + Y[2][2]*IX[0][2] + Y[3][2]*IX[0][3];

	M[1][0] = Y[0][0]*IX[1][0] + Y[1][0]*IX[1][1] + Y[2][0]*IX[1][2] + Y[3][0]*IX[1][3];
	M[1][1] = Y[0][1]*IX[1][0] + Y[1][1]*IX[1][1] + Y[2][1]*IX[1][2] + Y[3][1]*IX[1][3];
	M[1][2] = Y[0][2]*IX[1][0] + Y[1][2]*IX[1][1] + Y[2][2]*IX[1][2] + Y[3][2]*IX[1][3];

	M[2][0] = Y[0][0]*IX[2][0] + Y[1][0]*IX[2][1] + Y[2][0]*IX[2][2] + Y[3][0]*IX[2][3];
	M[2][1] = Y[0][1]*IX[2][0] + Y[1][1]*IX[2][1] + Y[2][1]*IX[2][2] + Y[3][1]*IX[2][3];
	M[2][2] = Y[0][2]*IX[2][0] + Y[1][2]*IX[2][1] + Y[2][2]*IX[2][2] + Y[3][2]*IX[2][3];

	M[3][0] = Y[0][0]*ix30     + Y[1][0]*ix31     + Y[2][0]*ix32     + Y[3][0]*ix33;
	M[3][1] = Y[0][1]*ix30     + Y[1][1]*ix31     + Y[2][1]*ix32     + Y[3][1]*ix33;
	M[3][2] = Y[0][2]*ix30     + Y[1][2]*ix31     + Y[2][2]*ix32     + Y[3][2]*ix33;

	*pix30 = ix30;
}

static void invertM(REAL M[4][3], REAL IM[4][3], REAL *d)
{
	REAL id;
	*d = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])+
	     M[0][1]*(M[1][2]*M[2][0]-M[1][0]*M[2][2])+
	     M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);

	id = 1.0/(*d);
	IM[0][0] = (M[1][1]*M[2][2]-M[1][2]*M[2][1])*id;
	IM[0][1] = (M[0][2]*M[2][1]-M[0][1]*M[2][2])*id;
	IM[0][2] = (M[0][1]*M[1][2]-M[0][2]*M[1][1])*id;

	IM[1][0] = (M[1][2]*M[2][0]-M[1][0]*M[2][2])*id;
	IM[1][1] = (M[0][0]*M[2][2]-M[0][2]*M[2][0])*id;
	IM[1][2] = (M[0][2]*M[1][0]-M[0][0]*M[1][2])*id;

	IM[2][0] = (M[1][0]*M[2][1]-M[1][1]*M[2][0])*id;
	IM[2][1] = (M[0][1]*M[2][0]-M[0][0]*M[2][1])*id;
	IM[2][2] = (M[0][0]*M[1][1]-M[0][1]*M[1][0])*id;

	IM[3][0] = (M[1][0]*(M[3][1]*M[2][2]-M[2][1]*M[3][2])+
		    M[1][1]*(M[2][0]*M[3][2]-M[3][0]*M[2][2])+
		    M[1][2]*(M[3][0]*M[2][1]-M[2][0]*M[3][1]))*id;
	IM[3][1] = (M[0][0]*(M[2][1]*M[3][2]-M[3][1]*M[2][2])+
		    M[0][1]*(M[3][0]*M[2][2]-M[2][0]*M[3][2])+
		    M[0][2]*(M[2][0]*M[3][1]-M[3][0]*M[2][1]))*id;
	IM[3][2] = (M[0][0]*(M[3][1]*M[1][2]-M[1][1]*M[3][2])+
		    M[0][1]*(M[1][0]*M[3][2]-M[3][0]*M[1][2])+
		    M[0][2]*(M[3][0]*M[1][1]-M[1][0]*M[3][1]))*id;
}

static void scan_line(REAL lin[2], int y, int z, int *n, int vox[][3], int maxv)
{
	REAL p[2], t;
	int x, xe;

	/* sort p into ascending order of x */
	p[0] = lin[0]; p[1] = lin[1];
	if (p[1]<p[0]) {t = p[1]; p[1] = p[0]; p[0] = t;}

	/* find voxels where x is integer */
	for(x=ceil(p[0]), xe=floor(p[1]); x<=xe; x++)
	{
		if ((*n)>=maxv-1)
		{
			(void)exit(1);
			mexErrMsgTxt("Too many voxels inside a tetrahedron");
		}
		vox[*n][0] = x;
		vox[*n][1] = y;
		vox[*n][2] = z;
		(*n)++;
	}
}

static void scan_triangle(REAL tri[][2], int z, int *n, int vox[][3], int maxv)
{
	REAL *p[3], *t, lin[2];
	REAL x1, x2, y1, y2;
	int y, ye, i;

	/* sort p into ascending order of y */
	p[0] = tri[0]; p[1] = tri[1]; p[2] = tri[2];
	if (p[1][1]<p[0][1]) {t = p[1]; p[1] = p[0]; p[0] = t;}
	if (p[2][1]<p[1][1]) {t = p[2]; p[2] = p[1]; p[1] = t;}
	if (p[1][1]<p[0][1]) {t = p[1]; p[1] = p[0]; p[0] = t;}

	/* find lower lines cutting triangle where y is integer */
	for(y=ceil(p[0][1]), ye=floor(p[1][1]); y<=ye; y++)
	{
		x1 = p[0][0]; y1 = p[0][1];
		for(i=0; i<2; i++)
		{
			x2 = p[i+1][0]; y2 = p[i+1][1];
			if (y2-y1<=0)
				lin[i] = (x1+x2)/2.0;
			else
				lin[i] = (x1*(y2-y)+x2*(y-y1))/(y2-y1);
		}
		scan_line(lin,y,z, n,vox,maxv);
	}

	/* find upper lines cutting triangle where y is integer */
	for(y=ceil(p[1][1]), ye=floor(p[2][1]); y<=ye; y++)
	{
		x2 = p[2][0]; y2 = p[2][1];
		for(i=0; i<2; i++)
		{
			x1 = p[i][0]; y1 = p[i][1];
			if (y2-y1<=0)
				lin[i] = (x1+x2)/2.0;
			else
				lin[i] = (x1*(y2-y)+x2*(y-y1))/(y2-y1);
		}
		scan_line(lin,y,z, n,vox,maxv);
	}
}


static void scan_tetrahedron(REAL Y[4][3], int *n, int vox[][3], int maxv)
{
	REAL *p[4], *t;
	REAL tri[4][2];
	REAL x1, x2, y1, y2, z1, z2;
	int z, ze, i;

	*n = 0;

	/* sort p into ascending order of z */
	p[0] = Y[0]; p[1] = Y[1]; p[2] = Y[2]; p[3] = Y[3];
	if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}
	if (p[2][2]<p[1][2]) {t = p[2]; p[2] = p[1]; p[1] = t;}
	if (p[3][2]<p[2][2]) {t = p[3]; p[3] = p[2]; p[2] = t;}
	if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}
	if (p[2][2]<p[1][2]) {t = p[2]; p[2] = p[1]; p[1] = t;}
	if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}

	/* find lower triangles that intersect tetrahedron where z is integer */
	for(z=ceil(p[0][2]), ze=floor(p[1][2]); z<=ze; z++)
	{
		x1 = p[0][0]; y1 = p[0][1]; z1 = p[0][2];
		for(i=0; i<3; i++)
		{
			x2 = p[i+1][0]; y2 = p[i+1][1]; z2 = p[i+1][2];
			if (z2-z1<=0)
			{
				tri[i][0] = (x1+x2)/2.0;
				tri[i][1] = (y1+y2)/2.0;
			}
			else
			{
				REAL t2 = z2-z, t1 = z-z1, t = z2-z1;
				tri[i][0] = (x1*t2+x2*t1)/t;
				tri[i][1] = (y1*t2+y2*t1)/t;
			}
		}
		scan_triangle(tri,z, n,vox,maxv);
	}

	/* find quadrilaterals that intersect tetrahedron where z is integer */
	/* each quadrilateral divided into two triangles */
	for(z=ceil(p[1][2]), ze=floor(p[2][2]); z<=ze; z++)
	{
		static int ii[] = {0,1,1,0}, jj[] = {3,3,2,2};

		for(i=0; i<4; i++)
		{
			x1 = p[ii[i]][0]; y1 = p[ii[i]][1]; z1 = p[ii[i]][2];
			x2 = p[jj[i]][0]; y2 = p[jj[i]][1]; z2 = p[jj[i]][2];
			if (z2-z1<=0)
			{
				tri[i][0] = (x1+x2)/2.0;
				tri[i][1] = (y1+y2)/2.0;
			}
			else
			{
				REAL t2 = z2-z, t1 = z-z1, t = z2-z1;
				tri[i][0] = (x1*t2+x2*t1)/t;
				tri[i][1] = (y1*t2+y2*t1)/t;
			}
		}
		scan_triangle(tri,z, n,vox,maxv);
		tri[1][0] = tri[3][0];
		tri[1][1] = tri[3][1];
		scan_triangle(tri,z, n,vox,maxv);
	}

	/* find upper triangles that intersect tetrahedron where z is integer */
	for(z=ceil(p[2][2]), ze=floor(p[3][2]); z<=ze; z++)
	{
		x2 = p[3][0]; y2 = p[3][1]; z2 = p[3][2];
		for(i=0; i<3; i++)
		{
			x1 = p[i][0]; y1 = p[i][1]; z1 = p[i][2];
			if (z2-z1<=0)
			{
				tri[i][0] = (x1+x2)/2.0;
				tri[i][1] = (y1+y2)/2.0;
			}
			else
			{
				REAL t2 = z2-z, t1 = z-z1, t = z2-z1;
				tri[i][0] = (x1*t2+x2*t1)/t;
				tri[i][1] = (y1*t2+y2*t1)/t;
			}
		}
		scan_triangle(tri,z, n,vox,maxv);
	}
}

static REAL x[2][5][4][3] = {
{{{ 0,0,0},{ 1,0,1},{ 1,0,0},{ 1,1,0}},
 {{ 0,0,0},{ 1,0,1},{ 0,1,1},{ 0,0,1}},
 {{ 0,0,0},{ 0,1,0},{ 0,1,1},{ 1,1,0}},
 {{ 0,0,0},{ 1,0,1},{ 1,1,0},{ 0,1,1}},
 {{ 1,1,1},{ 1,1,0},{ 0,1,1},{ 1,0,1}},},

{{{ 1,0,0},{ 0,0,1},{ 0,0,0},{ 0,1,0}},
 {{ 1,0,0},{ 0,0,1},{ 1,1,1},{ 1,0,1}},
 {{ 1,0,0},{ 1,1,0},{ 1,1,1},{ 0,1,0}},
 {{ 1,0,0},{ 0,0,1},{ 0,1,0},{ 1,1,1}},
 {{ 0,1,1},{ 0,1,0},{ 1,1,1},{ 0,0,1}},},
};

static REAL ix[2][5][4][4], dix[2][5][4][3];
static int off[2][4][5];

static void setup_consts(int dim[3])
{
	int i, k;
	for(k=0; k<2; k++)
		for(i=0; i<5; i++)
		{
			REAL dtx;
			int j;
			invertX(x[k][i], ix[k][i], dix[k][i], &dtx);
			for(j=0; j<4; j++)
				off[k][j][i] = x[k][i][j][0]+dim[0]*(x[k][i][j][1]+dim[1]*x[k][i][j][2]);
		}
}

static void invert_it(int x0, int x1, int x2, float *y0, float *y1, float *y2, int dim_f[3], float *iy0, float *iy1, float *iy2)
{
	int i, k;
	REAL d, ix30;
	k = (x0 + x1 + x2)%2;
	for(i=0; i<5; i++)
	{
		REAL Y[4][3], M[4][3], IM[4][3];
		int j, vox[MAXV][3], nvox;

		Y[0][0] = y0[off[k][0][i]]; Y[0][1] = y1[off[k][0][i]]; Y[0][2] = y2[off[k][0][i]];
		Y[1][0] = y0[off[k][1][i]]; Y[1][1] = y1[off[k][1][i]]; Y[1][2] = y2[off[k][1][i]];
		Y[2][0] = y0[off[k][2][i]]; Y[2][1] = y1[off[k][2][i]]; Y[2][2] = y2[off[k][2][i]];
		Y[3][0] = y0[off[k][3][i]]; Y[3][1] = y1[off[k][3][i]]; Y[3][2] = y2[off[k][3][i]];

		getM(Y, ix[k][i], dix[k][i], M, x0, x1, x2, &ix30);
		invertM(M, IM, &d);

		scan_tetrahedron(Y, &nvox, vox, MAXV);
		for(j=0; j<nvox; j++)
		{
/* printf("%d %d %d\n", vox[j][0], vox[j][1], vox[j][2]); */

			if ((vox[j][0]>=1) && (vox[j][0]<=dim_f[0]) &&
			    (vox[j][1]>=1) && (vox[j][1]<=dim_f[1]) &&
			    (vox[j][2]>=1) && (vox[j][2]<=dim_f[2]))
			{
				int o  = vox[j][0]+dim_f[0]*(vox[j][1]+dim_f[1]*vox[j][2]);
				iy0[o] = IM[0][0]*vox[j][0] + IM[1][0]*vox[j][1] + IM[2][0]*vox[j][2] + IM[3][0];
				iy1[o] = IM[0][1]*vox[j][0] + IM[1][1]*vox[j][1] + IM[2][1]*vox[j][2] + IM[3][1];
				iy2[o] = IM[0][2]*vox[j][0] + IM[1][2]*vox[j][1] + IM[2][2]*vox[j][2] + IM[3][2];
printf("%d %d %d\n", vox[j][0], vox[j][1], vox[j][2]);
			}
		}
	}
}



static void invert_field(int dim_g[3], float  y0[], float  y1[], float  y2[],
		         int dim_f[3], float iy0[], float iy1[], float iy2[])
{
	int x2, x1, x0;
	setup_consts(dim_g);

	y0  -= 1+dim_g[0]*(1 + dim_g[1]);
	y1  -= 1+dim_g[0]*(1 + dim_g[1]);
	y2  -= 1+dim_g[0]*(1 + dim_g[1]);

	iy0 -= 1+dim_f[0]*(1 + dim_f[1]);
	iy1 -= 1+dim_f[0]*(1 + dim_f[1]);
	iy2 -= 1+dim_f[0]*(1 + dim_f[1]);

	for(x2=1; x2<dim_g[2]; x2++)
	{
		for(x1=1; x1<dim_g[1]; x1++)
			for(x0=1; x0<dim_g[0]; x0++)
			{
				int o = x0 + dim_g[0]*(x1 + x2*dim_g[1]);
				invert_it(x0, x1, x2, y0+o, y1+o, y2+o, dim_f, iy0, iy1, iy2);
			}
		/* (void)printf(".");(void)fflush(stdout); */
	}
	/* (void)printf("\n");(void)fflush(stdout); */
}

static void setnan(float *dat, int n)
{
	int j;
	float NaN;
	NaN = mxGetNaN();
	for (j=0; j<n; j++)
		dat[j] = NaN;
}


static float *get_volume(const mxArray *ptr, int dims[3])
{
	int nd, i;
	const int *ldims;
	if (mxIsStruct(ptr) || !mxIsNumeric(ptr) || mxIsComplex(ptr) ||
                mxIsSparse(ptr) || !mxIsSingle(ptr))
		mexErrMsgTxt("Data must be a single precision floating point multi-dimensional array.");

	nd = mxGetNumberOfDimensions(ptr);
	if (nd>3)
		mexErrMsgTxt("Too many dimensions in data.");

	ldims = mxGetDimensions(ptr);
	for(i=0; i<nd; i++)
		dims[i] = ldims[i];
	for(i=nd; i<3; i++)
		dims[i] = 1;

	return((float *)mxGetPr(ptr));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	float  *y0=0,  *y1=0,  *y2=0, *iy0=0, *iy1=0, *iy2=0;
	int dim_g[3], dim_f[3];

        if (nrhs != 4 || nlhs > 3)
                mexErrMsgTxt("Inappropriate usage.");

	y0 = get_volume(prhs[0], dim_g);
	y1 = get_volume(prhs[1], dim_f);
	if (dim_g[0] != dim_f[0] || dim_g[1] != dim_f[1] || dim_g[2] != dim_f[2])
		mexErrMsgTxt("Incompatible dimensions.");
	y2 = get_volume(prhs[2], dim_f);
	if (dim_g[0] != dim_f[0] || dim_g[1] != dim_f[1] || dim_g[2] != dim_f[2])
		mexErrMsgTxt("Incompatible dimensions.");

	if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) ||
		mxIsComplex(prhs[3]) || !mxIsDouble(prhs[3]) || mxGetM(prhs[3]) * mxGetN(prhs[3]) != 3)
		mexErrMsgTxt("Output dimensions must be numeric, real, full, double and contain three elements.");

	dim_f[0] = mxGetPr(prhs[3])[0];
	dim_f[1] = mxGetPr(prhs[3])[1];
	dim_f[2] = mxGetPr(prhs[3])[2];

	plhs[0] = mxCreateNumericArray(3, dim_f,mxSINGLE_CLASS,mxREAL);
	plhs[1] = mxCreateNumericArray(3, dim_f,mxSINGLE_CLASS,mxREAL);
	plhs[2] = mxCreateNumericArray(3, dim_f,mxSINGLE_CLASS,mxREAL);

	iy0 = (float *)mxGetPr(plhs[0]);
	iy1 = (float *)mxGetPr(plhs[1]);
	iy2 = (float *)mxGetPr(plhs[2]);

	setnan(iy0, dim_f[0]*dim_f[1]*dim_f[2]);
	setnan(iy1, dim_f[0]*dim_f[1]*dim_f[2]);
	setnan(iy2, dim_f[0]*dim_f[1]*dim_f[2]);

	invert_field(dim_g, y0, y1, y2, dim_f, iy0, iy1, iy2);
}
