#ifndef lint
static char sccsid[]="%W% John Ashburner FIL %E%";
#endif

/* Do a silent deletion of files on disk */

#include "cmex.h"
#include <unistd.h>

#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	int i;
	if (nlhs != 0) mexErrMsgTxt("Too many output arguments.");

	for(i=0; i<nrhs; i++)
	{
		Matrix *matptr;
		matptr = prhs[i];
		if (mxIsString(matptr))
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
