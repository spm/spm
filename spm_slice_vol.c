#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
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
	int m,n, k, hold, xdim, ydim, zdim, status;
	double *mat, *ptr, *img, background=0.0;

	if (nrhs != 4 || nlhs > 1)
	{
		mexErrMsgTxt("Inappropriate usage.");
	}

	map = get_map(prhs[0]);

	xdim = abs((int)map->xdim);
	ydim = abs((int)map->ydim);
	zdim = abs((int)map->zdim);

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
	m = abs((int)ptr[0]);
	n = abs((int)ptr[1]);
	plhs[0] = mxCreateFull(m,n,REAL);
	img = mxGetPr(plhs[0]);

	if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 1 && mxGetM(prhs[3])*mxGetN(prhs[3]) != 2)
	{
		mexErrMsgTxt("Hold & background argument must have one or two element(s).");
	}
	hold = abs((int)(*(mxGetPr(prhs[3]))));
	if (hold > 127 || hold == 2)
		mexErrMsgTxt("Bad hold value.");

	if (mxGetM(prhs[3])*mxGetN(prhs[3]) > 1)
		background = mxGetPr(prhs[3])[1];

	status = 1;
	if (map->datatype == UNSIGNED_CHAR)
		status = slice_uchar(mat, img, m, n, map->data, xdim, ydim, zdim,
			hold, background, map->scalefactor);
	else if (map->datatype == SIGNED_SHORT)
		status = slice_short(mat, img, m, n, map->data, xdim, ydim, zdim,
			hold, background, map->scalefactor);
	else if (map->datatype == SIGNED_INT)
		status = slice_int(mat, img, m, n, map->data, xdim, ydim, zdim,
			hold, background, map->scalefactor);
	else if (map->datatype == FLOAT)
		status = slice_float(mat, img, m, n, map->data, xdim, ydim, zdim,
			hold, background, map->scalefactor);
	else if (map->datatype == DOUBLE)
		status = slice_double(mat, img, m, n, map->data, xdim, ydim, zdim,
			hold, background, map->scalefactor);
	if (status)
		mexErrMsgTxt("Slicing failed.");
}
