#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
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
	int m,n, k, hold, xdim,ydim,zdim;
	double background=0.0;

	if (nrhs != 5 || nlhs > 4)
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
		(mxGetM(prhs[4])*mxGetN(prhs[4]) != 1 && mxGetM(prhs[4])*mxGetN(prhs[4]) != 2))
		mexErrMsgTxt("Bad hold & background argument.");

	hold = (int)(*(mxGetPr(prhs[4])));

	if (abs(hold) > 127)
		mexErrMsgTxt("Bad hold value.");

	if (mxGetM(prhs[4])*mxGetN(prhs[4]) > 1)
		background = mxGetPr(prhs[4])[1];

	if (nlhs<=1)
	{
		plhs[0] = mxCreateFull(m,n,REAL);

		resample(m*n, map->data, mxGetPr(plhs[0]),
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]), xdim, ydim, zdim,
			hold, background, map->scalefactor,0.0,map->datatype);
	}
	else
	{
		if (hold==0)
			mexErrMsgTxt("This wont work for nearest neighbour resampling.");
		plhs[0] = mxCreateFull(m,n,REAL);
		plhs[1] = mxCreateFull(m,n,REAL);
		plhs[2] = mxCreateFull(m,n,REAL);
		plhs[3] = mxCreateFull(m,n,REAL);

		resample_d(m*n, map->data, mxGetPr(plhs[0]),mxGetPr(plhs[1]),mxGetPr(plhs[2]),mxGetPr(plhs[3]),
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]), xdim, ydim, zdim,
			hold, background, map->scalefactor,0.0,map->datatype);
	}
}
