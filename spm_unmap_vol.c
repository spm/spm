#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include "mex.h"
#include "spm_map.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray *matrix_ptr;
	MAP *current_map;
	int j, n;

	if (nrhs != 1 || nlhs > 0)
		mexErrMsgTxt("Inappropriate usage.");

	matrix_ptr = (mxArray *)prhs[0];
	if ((!mxIsNumeric(matrix_ptr) || mxIsComplex(matrix_ptr) ||
		mxIsSparse(matrix_ptr) || !mxIsDouble(matrix_ptr) ||
		mxGetM(matrix_ptr) != (sizeof(MAP)+sizeof(double)-1)/sizeof(double)))

	n = mxGetN(matrix_ptr);

	for (j=0; j<n; j++)
	{
		current_map = (MAP *)(mxGetPr(matrix_ptr)
			+ j*((sizeof(MAP)+sizeof(double)-1)/sizeof(double)));
		if (current_map->magic != MAGIC)
			mexErrMsgTxt("Bad magic number in image handle.");
		if (current_map->pid != getpid())
			mexErrMsgTxt("Invalid image handle (from old session).");

		(void)munmap(current_map->map, current_map->off+current_map->len);
		current_map->map=0;
		current_map->off=0;
		current_map->len=0;
		current_map->magic=0;
		current_map->pid=0;
	}
}
