#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include <sys/types.h>
#include <unistd.h>
#ifdef SPM_WIN32
#include <process.h>
#endif
#include "spm_vol_access.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MAP *current_map;
	int j, n;

	if (nrhs != 1 || nlhs > 0)
		mexErrMsgTxt("Inappropriate usage.");

	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ||
		mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) ||
		mxGetM(prhs[0]) != (sizeof(MAP)+sizeof(double)-1)/sizeof(double))
		mexErrMsgTxt("Wrong sized matrix/vector.");

	n = mxGetN(prhs[0]);

	for (j=0; j<n; j++)
	{
		current_map = (MAP *)mxGetPr(prhs[0]) +j;
		if (current_map->magic != MAGIC)
			mexErrMsgTxt("Bad magic number in image handle.");
		if (current_map->pid != getpid())
			mexErrMsgTxt("Invalid image handle (from old session).");
		
		#ifdef SPM_WIN32
		(void)unmap_file(current_map->map);
		#else
		(void)munmap(current_map->map, current_map->off+current_map->len);
		#endif

		current_map->map=0;
		current_map->off=0;
		current_map->len=0;
		current_map->magic=0;
		current_map->pid=0;
	}
}
