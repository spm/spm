#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include <math.h>
#include "mex.h"
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
	hold = (int)(*(mxGetPr(prhs[3])));
	if (abs(hold) > 127)
		mexErrMsgTxt("Bad hold value.");

	if (mxGetM(prhs[3])*mxGetN(prhs[3]) > 1)
		background = mxGetPr(prhs[3])[1];

	status = slice(mat, img, m, n, map->data, xdim, ydim, zdim,
		hold, background, map->scalefactor,0.0,map->datatype);
	if (status)
		mexErrMsgTxt("Slicing failed.");
}
