#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif
/*
Extracts x, y and z co-ordinates from an image
*/

#include "spm_mapping.h"

static void transverse_plane(double img[], MAPTYPE *map, int i)
{
	static double mat[16] = {
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0};
	mat[14] = i+1.0;
	(void)slice(mat, img, map->dim[0], map->dim[1], map, 0, 0);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MAPTYPE *map, *get_maps();
	int m,n,k,i, nxyz;
	double *XYZ, *img;
	unsigned char *vol;

	if (nrhs != 1 || nlhs > 1)
	{
		mexErrMsgTxt("Inappropriate usage.");
	}

	map = get_maps(prhs[0], &n);
	if (n!=1)
	{
		free_maps(map, n);
		mexErrMsgTxt("Bad image handle dimensions.");
	}

	m = map->dim[0];
	n = map->dim[1];
	k = map->dim[2];

	vol  = (unsigned char *)mxCalloc(m*n*k,sizeof(unsigned char));

	img  = (double *)mxCalloc(m*n,sizeof(double));
	for(i=0, nxyz=0; i<k; i++)
	{
		double *ptr, *ptend;
		unsigned char *vp;
		transverse_plane(img,map,i);
		
		for(ptr=img, ptend=img+m*n, vp=vol+i*m*n; ptr<ptend; ptr++, vp++)
			if (mxIsFinite(*ptr) && *ptr)
			{
				*vp = 1;
				nxyz++;
			}
	}
	(void)mxFree((char *)img);
	free_maps(map, 1);

	plhs[0] = mxCreateDoubleMatrix(3,nxyz,mxREAL);
	XYZ     = mxGetPr(plhs[0]);

	for(i=0; i<k; i++)
	{
		int i1, j1;
		unsigned char *vp;

		for(j1=0, vp=vol+i*m*n; j1<n; j1++)
			for(i1=0; i1<m; i1++, vp++)
				if (*vp)
				{
					/* voxel coordinates...
					XYZ[0] = i1+1.0;
					XYZ[1] = j1+1.0;
					XYZ[2] = i +1.0;
					*/

					XYZ[0] = map->mat[0+0*4]*(i1+1) + map->mat[0+1*4]*(j1+1) +
					         map->mat[0+2*4]*( i+1) + map->mat[0+3*4];
					XYZ[1] = map->mat[1+0*4]*(i1+1) + map->mat[1+1*4]*(j1+1) +
					         map->mat[1+2*4]*( i+1) + map->mat[1+3*4];
					XYZ[2] = map->mat[2+0*4]*(i1+1) + map->mat[2+1*4]*(j1+1) +
					         map->mat[2+2*4]*( i+1) + map->mat[2+3*4];
					XYZ+=3;
				}
	}
	(void)mxFree((char *)vol);
}
