#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif

#include <math.h>
#include "cmex.h"
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
	int m,n, k, hold, xdim,ydim,zdim;
	double *img;

	if (nrhs != 5 || nlhs > 1)
	{
		mexErrMsgTxt("Inappropriate usage.");
	}

	map=get_map(prhs[0]);

	xdim = abs((int)map->xdim);
	ydim = abs((int)map->ydim);
	zdim = abs((int)map->zdim);

	for(k=1; k<=3; k++)
		if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
			!mxIsFull(prhs[k]) || !mxIsDouble(prhs[k]))
			mexErrMsgTxt("Coordinates must be numeric, real, full and double.");

	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if (mxGetM(prhs[2]) != m || mxGetN(prhs[2]) != n ||
		mxGetM(prhs[3]) != m || mxGetN(prhs[3]) != n)
		mexErrMsgTxt("Coordinates must have compatible dimensions.");

	if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) ||
		!mxIsFull(prhs[4]) || !mxIsDouble(prhs[4]) ||
		mxGetM(prhs[4])*mxGetN(prhs[4]) != 1)
		mexErrMsgTxt("Bad hold argument.");

	hold = abs((int)(*(mxGetPr(prhs[4]))));

	plhs[0] = mxCreateFull(m,n,REAL);
	img = mxGetPr(plhs[0]);

	if (hold > 127 || hold == 2)
		mexErrMsgTxt("Bad hold value.");

	if (map->datatype == UNSIGNED_CHAR)
		resample_uchar(m*n, map->data, img,
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]), xdim, ydim, zdim, hold);
	else if (map->datatype == SIGNED_SHORT)
		resample_short(m*n, map->data, img,
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]), xdim, ydim, zdim, hold);
	else if (map->datatype == SIGNED_INT)
		resample_int(m*n, map->data, img,
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]), xdim, ydim, zdim, hold);
	else if (map->datatype == FLOAT)
		resample_float(m*n, map->data, img,
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]), xdim, ydim, zdim, hold);
	else if (map->datatype == DOUBLE)
		resample_double(m*n, map->data, img,
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]), xdim, ydim, zdim, hold);
	else
		mexErrMsgTxt("Bad datatype.");

	/* final rescale */
	if (map->scalefactor != 1.0 && map->scalefactor != 0)
		for (k=0; k<m*n; k++)
			img[k] *= map->scalefactor;
}
