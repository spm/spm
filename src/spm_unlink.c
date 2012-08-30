/*
 * $Id: spm_unlink.c 4872 2012-08-30 15:29:15Z guillaume $
 * John Ashburner
 */

/* Do a silent deletion of files on disk */

#include <stdio.h>
#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i;
    if (nlhs != 0) mexErrMsgTxt("Too many output arguments.");

    for(i=0; i<nrhs; i++)
    {
        const mxArray *matptr = prhs[i];
        if (mxIsChar(matptr))
        {
            char *str = NULL;
            mwIndex k;
            mwSize len;
            
            str = mxArrayToString(matptr);
            len = strlen(str);

            /* delete white space */
            for (k=len-1;k>=0;k--)
                if (str[k] == ' ')
                    str[k] = '\0';
                else
                    break;
            
            remove(str); /* not bothered about return status */
            mxFree(str);
        }
        else
            mexErrMsgTxt("Filename should be a string.");

    }
}
