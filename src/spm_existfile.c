/*
 * $Id: spm_existfile.c 4901 2012-09-05 15:10:48Z guillaume $
 * Guillaume Flandin
 */
 
#include "io64.h"
#include "mex.h"

#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int status     = 0;
    char *filename = NULL;
    structStat stbuf;
    
    if (nrhs != 1)
    {
        mexErrMsgTxt("One input only required.");
    }
    else
    {
        if (!mxIsChar(prhs[0]))
        {
            mexErrMsgTxt("Input must be a string.");
        }
        filename = mxArrayToString(prhs[0]);
        
        if ((getFileStat(filename, &stbuf) == 0) && (S_ISREG(stbuf.st_mode)))
        {
        	status = 1;
        }

        mxFree(filename);
    }
        
    plhs[0] = mxCreateLogicalScalar(status);
}
