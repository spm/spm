#include "spm_map.h"
#include "mex.h"
#include "spm_vol_utils.h"
#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

int get_datasize(int type)
{
	if (type == UNSIGNED_CHAR) return(8);
	if (type == SIGNED_SHORT || type == SIGNED_SHORT_S) return(16);
	if (type == SIGNED_INT || type == SIGNED_INT_S) return(32);
	if (type == FLOAT || type == FLOAT_S) return(32);
	if (type == DOUBLE || type == DOUBLE_S) return(64);
	return(0);
}

MAPTYPE *get_maps(mxArray *matrix_ptr, int *pn)
{
	MAP *current_map;
	MAPTYPE *maps;
	int n, j;

	if ((!mxIsNumeric(matrix_ptr) || mxIsComplex(matrix_ptr) ||
		mxIsSparse(matrix_ptr) || !mxIsDouble(matrix_ptr) ||
		mxGetM(matrix_ptr) != (sizeof(MAP)+sizeof(double)-1)/sizeof(double)))
	{
		mexErrMsgTxt("Bad image handle dimensions.");
	}
	n = mxGetN(matrix_ptr);

	for (j=0; j<n; j++)
	{
		int xdim, ydim, zdim, datasize;
		current_map = (MAP *)mxGetPr(matrix_ptr)
			+ j*((sizeof(MAP)+sizeof(double)-1)/sizeof(double));
		if (current_map->magic != MAGIC)
		{
			mexErrMsgTxt("Bad magic number in image handle.");
		}
		if (current_map->pid != getpid())
		{
			mexErrMsgTxt("Invalid image handle (from old session).");
		}
		datasize = get_datasize(current_map->dtype);	
		if (datasize == 0)
			mexErrMsgTxt("Bad datatype in image handle.");

		xdim = abs((int)current_map->dim[0]);
		ydim = abs((int)current_map->dim[1]);
		zdim = abs((int)current_map->dim[2]);
		if ((xdim*ydim*zdim*datasize+7)/8 != current_map->len)
		{
			mexErrMsgTxt("Who screwed up the handle??");
		}
	}

	maps = (MAPTYPE *)mxCalloc(n, sizeof(MAPTYPE));
	for (j=0; j<n; j++)
	{
		current_map = (MAP *)(mxGetPr(matrix_ptr)
			+ j*((sizeof(MAP)+sizeof(double)-1)/sizeof(double)));

		maps[j].dim[0]    = current_map->dim[0];
		maps[j].dim[1]    = current_map->dim[1];
		maps[j].dim[2]    = current_map->dim[2];

		maps[j].pixdim[0] = current_map->pixdim[0];
		maps[j].pixdim[1] = current_map->pixdim[1];
		maps[j].pixdim[2] = current_map->pixdim[2];

		maps[j].scale     = current_map->scale;
		maps[j].offset    = 0.0;

		maps[j].dtype     = current_map->dtype;
		maps[j].data      = current_map->data;
	}
	*pn = n;
	return(maps);
}

void free_maps(MAPTYPE *maps)
{
	mxFree((char *)maps);
}
