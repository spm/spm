/*
 * $Id: spm_unlink.c 112 2005-05-04 18:20:52Z john $
 */

/* Do a silent deletion of files on disk */

#include "mex.h"
#ifndef SPM_WIN32
#include <unistd.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i;
	if (nlhs != 0) mexErrMsgTxt("Too many output arguments.");

	for(i=0; i<nrhs; i++)
	{
		const mxArray *matptr;
		matptr = prhs[i];
		if (!mxIsNumeric(matptr))
		{
			char *str;
			int k, stlen;
			
			stlen = mxGetN(matptr);
			str = (char *)mxCalloc(stlen+1, sizeof(char));
			mxGetString(matptr,str,stlen+1);

			/* delete white space */
			for(k=0; k<stlen; k++)
				if (str[k] == ' ')
				{
					str[k] = '\0';
					break;
				}
			unlink(str); /* not bothered about return status */
			mxFree(str);
		}
		else
			mexErrMsgTxt("Filename should be a string.");

	}
}
