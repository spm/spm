#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

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
	int i, n;
	double *dat;
	if (nrhs != 1 || nlhs > 0)
		mexErrMsgTxt("Inappropriate usage.");

	map=get_map(prhs[0]);
	(void)munmap(map->map, map->off+map->len);

	dat = mxGetPr(prhs[0]);
	n = mxGetN(prhs[0])*mxGetM(prhs[0]);
	for(i=0; i< n; i++)
		dat[i] = 0.0;
}
