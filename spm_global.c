/* returns the global mean for a memory mapped volume image
  FORMAT [G] = spm_global(V)
  V   - memory mapped volume
  G   - mean global activity
 ____________________________________________________________________________
 
  spm_global returns the mean counts integrated over all the  
  slices from the volume
 
  The mean is estimated after discounting voxels outside the object
  using a criteria of greater than > (global mean)/8
*/

#ifndef lint
static char sccsid[]="%W% anon %E%";
#endif

#include <math.h>
#include "cmex.h"
#include "volume.h"

#define macro\
	for(i=0, s1=0.0; i<m; i++)\
		s1=s1+d[i];\
	s1=s1/m/8.0;\
	for(i=0, s2=0.0, n=0.0; i<m; i++)\
		if (d[i]>s1){n++;s2+=d[i];}\
	s2=s2/n;


#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	int i, m, n;
	double s1=0.0, s2=0.0;
	MAPTYPE *map;

	if (nrhs != 1 || nlhs > 1)
	{
		mexErrMsgTxt("Inappropriate usage.");
	}

	map = get_map(prhs[0]);
	m = map->zdim*map->ydim*map->xdim;

	if (map->datatype == UNSIGNED_CHAR)
	{
		unsigned char *d = (unsigned char *)(map->data);
		macro
	}
	else if (map->datatype == SIGNED_SHORT)
	{
		short *d = (short *)(map->data);
		macro
	}
	else if (map->datatype == SIGNED_INT)
	{
		int *d = (int *)(map->data);
		macro
	}
	else if (map->datatype == FLOAT)
	{
		float *d = (float *)(map->data);
		macro
	}
	else if (map->datatype == DOUBLE)
	{
		double *d = (double *)(map->data);
		macro
	}
	plhs[0] = mxCreateFull(1,1,REAL);
	mxGetPr(plhs[0])[0]=s2*map->scalefactor;
}
